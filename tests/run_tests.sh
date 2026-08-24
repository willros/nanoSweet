#!/usr/bin/env bash
set -euo pipefail

PASS=0
FAIL=0
NANOMUX=${NANOMUX:-./nanomux}
TEST_TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEST_TMP_DIR"' EXIT

# ---------- helpers ----------

count_reads() {
    gunzip -c "$1" | grep -c '^@' || echo 0
}

assert_eq() {
    local desc="$1" expected="$2" actual="$3"
    if [ "$expected" = "$actual" ]; then
        PASS=$((PASS + 1))
    else
        FAIL=$((FAIL + 1))
        echo "  FAIL: $desc (expected '$expected', got '$actual')"
    fi
}

assert_file_exists() {
    local desc="$1" path="$2"
    if [ -f "$path" ]; then
        PASS=$((PASS + 1))
    else
        FAIL=$((FAIL + 1))
        echo "  FAIL: $desc — file not found: $path"
    fi
}

assert_file_not_exists() {
    local desc="$1" path="$2"
    if [ ! -f "$path" ]; then
        PASS=$((PASS + 1))
    else
        FAIL=$((FAIL + 1))
        echo "  FAIL: $desc — file should not exist: $path"
    fi
}

assert_read_in_output() {
    local desc="$1" read_name="$2" gz_file="$3"
    if gunzip -c "$gz_file" | grep -q "^@${read_name}$"; then
        PASS=$((PASS + 1))
    else
        FAIL=$((FAIL + 1))
        echo "  FAIL: $desc — read '$read_name' not found in $gz_file"
    fi
}

assert_read_not_in_output() {
    local desc="$1" read_name="$2" gz_file="$3"
    if ! gunzip -c "$gz_file" | grep -q "^@${read_name}$"; then
        PASS=$((PASS + 1))
    else
        FAIL=$((FAIL + 1))
        echo "  FAIL: $desc — read '$read_name' should not be in $gz_file"
    fi
}

get_match_count() {
    local csv="$1" barcode="$2"
    grep "^${barcode}," "$csv" | cut -d, -f2
}

# ---------- Test 1: Single barcode, k=0, no trim ----------
echo "TEST 1: Single barcode, k=0, no trim"
OUT="$TEST_TMP_DIR/test1"
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 0 -j 1 >/dev/null 2>&1

assert_eq "BC_A match count" "4" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_A")"
assert_eq "BC_B match count" "1" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_B")"

assert_read_in_output "read_a_fw_k0 in BC_A" "read_a_fw_k0" "$OUT/BC_A.fq.gz"
assert_read_in_output "read_a_3prime in BC_A" "read_a_3prime" "$OUT/BC_A.fq.gz"
assert_read_in_output "read_a_dual_fwd in BC_A" "read_a_dual_fwd" "$OUT/BC_A.fq.gz"
assert_read_in_output "read_a_dual_rev in BC_A" "read_a_dual_rev" "$OUT/BC_A.fq.gz"
assert_read_in_output "read_b_fw_k0 in BC_B" "read_b_fw_k0" "$OUT/BC_B.fq.gz"
assert_read_not_in_output "read_nomatch not in BC_A" "read_nomatch" "$OUT/BC_A.fq.gz"
assert_read_not_in_output "read_nomatch not in BC_B" "read_nomatch" "$OUT/BC_B.fq.gz"

# ---------- Test 2: Single barcode, k=1 ----------
echo "TEST 2: Single barcode, k=1"
OUT="$TEST_TMP_DIR/test2"
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 1 -j 1 >/dev/null 2>&1

assert_eq "BC_A match count with k=1" "5" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_A")"
assert_read_in_output "read_a_fw_k1 in BC_A with k=1" "read_a_fw_k1" "$OUT/BC_A.fq.gz"

# ---------- Test 3: Single barcode, k=0, with trim ----------
echo "TEST 3: Single barcode, k=0, with trim"
OUT="$TEST_TMP_DIR/test3"
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 0 -j 1 -t >/dev/null 2>&1

# Trimmed reads should be shorter than 200
trimmed_len=$(gunzip -c "$OUT/BC_A.fq.gz" | awk 'NR==2{print length($0)}')
assert_eq "trimmed read_a_fw_k0 length" "178" "$trimmed_len"

# 3' match: trimmed to before barcode position
trimmed_3prime=$(gunzip -c "$OUT/BC_A.fq.gz" | awk '/^@read_a_3prime/{getline; print length($0)}')
assert_eq "trimmed read_a_3prime length" "160" "$trimmed_3prime"

# ---------- Test 4: Dual barcode, k=0 ----------
echo "TEST 4: Dual barcode, k=0"
OUT="$TEST_TMP_DIR/test4"
$NANOMUX -b tests/test_barcodes_dual.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 0 -j 1 >/dev/null 2>&1

assert_eq "BC_A dual match count" "2" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_A")"
assert_eq "BC_B dual match count" "0" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_B")"
assert_read_in_output "read_a_dual_fwd in dual BC_A" "read_a_dual_fwd" "$OUT/BC_A.fq.gz"
assert_read_in_output "read_a_dual_rev in dual BC_A" "read_a_dual_rev" "$OUT/BC_A.fq.gz"

# Reads with only one barcode should NOT match in dual mode
assert_read_not_in_output "read_a_fw_k0 not in dual BC_A" "read_a_fw_k0" "$OUT/BC_A.fq.gz"

# ---------- Test 5: Multi-thread determinism (j=1 vs j=4) ----------
echo "TEST 5: Multi-thread determinism"
OUT1="$TEST_TMP_DIR/test5_j1"
OUT4="$TEST_TMP_DIR/test5_j4"
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT1" -p 50 -k 0 -j 1 >/dev/null 2>&1
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT4" -p 50 -k 0 -j 4 >/dev/null 2>&1

matches_j1=$(cat "$OUT1/nanomux_matches.csv")
matches_j4=$(cat "$OUT4/nanomux_matches.csv")
assert_eq "j=1 vs j=4 matches.csv identical" "$matches_j1" "$matches_j4"

# ---------- Test 6: Empty input ----------
echo "TEST 6: Empty input"
OUT="$TEST_TMP_DIR/test6"
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_empty.fastq -o "$OUT" -p 50 -k 0 -j 1 >/dev/null 2>&1

assert_eq "BC_A empty input" "0" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_A")"
assert_eq "BC_B empty input" "0" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_B")"
assert_file_not_exists "BC_A.fq.gz deleted when empty" "$OUT/BC_A.fq.gz"
assert_file_not_exists "BC_B.fq.gz deleted when empty" "$OUT/BC_B.fq.gz"

# ---------- Test 7: Short reads reported correctly ----------
echo "TEST 7: Short reads count"
OUT="$TEST_TMP_DIR/test7"
$NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 0 -j 1 2>&1 | grep -o "Reads shorter than p: [0-9]* reads" > "$TEST_TMP_DIR/short_msg.txt"
short_count=$(cat "$TEST_TMP_DIR/short_msg.txt" | grep -o '[0-9]*' | head -1)
assert_eq "reads shorter than p" "1" "$short_count"

# ---------- Test 8: Invalid k=4 ----------
echo "TEST 8: Invalid k=4 rejected"
OUT="$TEST_TMP_DIR/test8"
if $NANOMUX -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 4 -j 1 >/dev/null 2>&1; then
    FAIL=$((FAIL + 1))
    echo "  FAIL: k=4 should return non-zero exit code"
else
    PASS=$((PASS + 1))
fi

# ---------- Test 9: Reverse orientation in dual mode ----------
echo "TEST 9: Reverse orientation in dual mode"
OUT="$TEST_TMP_DIR/test9"
$NANOMUX -b tests/test_barcodes_dual.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 0 -j 1 >/dev/null 2>&1
assert_read_in_output "read_a_dual_rev matched in dual (rv...fw_comp)" "read_a_dual_rev" "$OUT/BC_A.fq.gz"

# ---------- Test 10: Output remains readable across multiple read buffers ----------
echo "TEST 10: Multi-buffer gzip output"
MULTI_BARCODES="$TEST_TMP_DIR/multi_buffer_barcodes.csv"
MULTI_FASTQ="$TEST_TMP_DIR/multi_buffer.fastq"
printf 'name,forward\nBC_MULTI,ACGT\n' > "$MULTI_BARCODES"
: > "$MULTI_FASTQ"
for i in $(seq 1 10005); do
    printf '@multi_%s\nACGTAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n' "$i" >> "$MULTI_FASTQ"
done

OUT="$TEST_TMP_DIR/test10"
$NANOMUX -b "$MULTI_BARCODES" -f "$MULTI_FASTQ" -o "$OUT" -p 12 -k 0 -j 2 >/dev/null 2>&1
assert_eq "multi-buffer match count" "10005" "$(get_match_count "$OUT/nanomux_matches.csv" "BC_MULTI")"
assert_eq "multi-buffer output read count" "10005" "$(count_reads "$OUT/BC_MULTI.fq.gz")"
if gzip -t "$OUT/BC_MULTI.fq.gz"; then
    PASS=$((PASS + 1))
else
    FAIL=$((FAIL + 1))
    echo "  FAIL: multi-buffer output is not a valid gzip stream"
fi

# ---------- Test 11: Many matching barcodes under a low descriptor limit ----------
echo "TEST 11: Low descriptor limit with many barcodes"
MANY_BARCODES="$TEST_TMP_DIR/many_barcodes.csv"
ONE_READ="$TEST_TMP_DIR/one_read.fastq"
printf 'name,forward\n' > "$MANY_BARCODES"
for i in $(seq 1 300); do
    printf 'BC_%03d,ACGT\n' "$i" >> "$MANY_BARCODES"
done
printf '@fd_read\nACGTAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$ONE_READ"

OUT="$TEST_TMP_DIR/test11"
if (ulimit -n 32; "$NANOMUX" -b "$MANY_BARCODES" -f "$ONE_READ" -o "$OUT" -p 12 -k 0 -j 64) >"$TEST_TMP_DIR/test11.stdout" 2>"$TEST_TMP_DIR/test11.stderr"; then
    PASS=$((PASS + 1))
else
    FAIL=$((FAIL + 1))
    echo "  FAIL: nanomux should succeed with 300 barcodes under ulimit -n 32"
    sed -n '1,20p' "$TEST_TMP_DIR/test11.stderr"
fi
if grep -q 'Concurrent gzip outputs limited to 16' "$TEST_TMP_DIR/test11.stderr" &&
   grep -q 'threads: 64' "$TEST_TMP_DIR/test11.stderr"; then
    PASS=$((PASS + 1))
else
    FAIL=$((FAIL + 1))
    echo "  FAIL: output handles were not capped while preserving all 64 workers"
fi
many_summary_rows=$(tail -n +2 "$OUT/nanomux_matches.csv" | wc -l | tr -d '[:space:]')
many_output_files=$(find "$OUT" -name 'BC_*.fq.gz' -type f | wc -l | tr -d '[:space:]')
assert_eq "all low-limit barcodes summarized" "300" "$many_summary_rows"
assert_eq "all low-limit barcode outputs created" "300" "$many_output_files"
assert_eq "low-limit output remains readable" "1" "$(count_reads "$OUT/BC_001.fq.gz")"

# ---------- Test 12: Summary open failure exits cleanly and flushes outputs ----------
echo "TEST 12: Summary open failure cleanup"
OUT="$TEST_TMP_DIR/test12"
if (ulimit -n 5; "$NANOMUX" -b tests/test_barcodes_single.csv -f tests/test_known.fastq -o "$OUT" -p 50 -k 0 -j 1) >"$TEST_TMP_DIR/test12.stdout" 2>"$TEST_TMP_DIR/test12.stderr"; then
    summary_failure_status=0
else
    summary_failure_status=$?
fi
assert_eq "summary open failure returns a normal error" "1" "$summary_failure_status"
assert_eq "output is flushed before summary failure exit" "4" "$(count_reads "$OUT/BC_A.fq.gz")"
if gzip -t "$OUT/BC_A.fq.gz"; then
    PASS=$((PASS + 1))
else
    FAIL=$((FAIL + 1))
    echo "  FAIL: output gzip is invalid after summary open failure"
fi
if grep -q 'Could NOT create summary file.*Too many open files' "$TEST_TMP_DIR/test12.stderr"; then
    PASS=$((PASS + 1))
else
    FAIL=$((FAIL + 1))
    echo "  FAIL: summary open failure did not report its path and system error"
fi
if grep -q 'Segmentation fault' "$TEST_TMP_DIR/test12.stderr"; then
    FAIL=$((FAIL + 1))
    echo "  FAIL: summary open failure caused a segmentation fault"
else
    PASS=$((PASS + 1))
fi

# ---------- Summary ----------
echo ""
echo "=== Integration Tests: $PASS passed, $FAIL failed ==="
exit $FAIL
