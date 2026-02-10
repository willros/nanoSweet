# nanomux speed-up benchmarks

Branch: `speed-up`
Date: 2026-02-10
System: Darwin 24.2.0 (Apple Silicon)
Compiler: `cc -O3`

## Test setup

- Input: `tests/test.fastq` concatenated 100x → 478,800 reads
- Barcodes: `tests/bc_test.csv` (9 dual barcodes)
- Barcode position: `-p 100`

## Results

| Scenario   | Old (main) | New (speed-up) | Speedup |
|------------|------------|----------------|---------|
| k=0, j=1  | 26.2s      | 4.6s           | 5.7x    |
| k=0, j=4  | 10.4s      | 3.5s           | 3.0x    |
| k=1, j=1  | 25.6s      | 22.9s          | 1.1x    |
| k=1, j=4  | 10.3s      | 8.4s           | 1.2x    |

## Changes

### 1. Levenshtein k=0 fast path (common.h)

When k=0 (exact matching, the most common case), skip the O(n*m) DP matrix
entirely and use a simple `memcmp` substring scan. This is the main source of
the ~6x speedup at k=0.

### 2. Row-major DP with early return (common.h)

For k>0, keep the original row-major loop order (good cache locality) but add
an inline early return on the last row: as soon as `dp[needle_len][j] <= k`,
return immediately instead of filling the rest of the matrix.

Note: column-major fill was tested first (early return after each column) but
was 2x slower due to cache stride — the inner loop accessed elements
`haystack_len+1` apart instead of contiguous memory. Reverted to row-major.

### 3. Read-parallel threading (nanomux.c)

**Before (barcode-parallel):** one thread per barcode, each scans all reads.
With 9 barcodes and 4 threads, 3 sequential rounds. Every read checked against
every barcode regardless of matches.

**After (read-parallel):** `num_threads` threads, each processes a chunk of
reads against all barcodes. Early exit on first match — once a read matches a
barcode, remaining barcodes are skipped.

This also fixes a correctness issue: previously a read could be written to
multiple barcode output files. Now each read is assigned to at most one barcode
(first match in CSV order).

Requires one `pthread_mutex_t` per barcode to protect concurrent gz writes.
Contention is low since matches are sparse.
