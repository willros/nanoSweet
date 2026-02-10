#define COMMON_IMPLEMENTATION
#include "../common.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (cond) { \
        tests_passed++; \
    } else { \
        tests_failed++; \
        printf("  FAIL: %s (line %d)\n", msg, __LINE__); \
    } \
} while(0)

#define TEST(name) printf("TEST: %s\n", name)

// ---- complement ----
void test_complement(void) {
    TEST("complement");
    ASSERT(complement('A') == 'T', "A -> T");
    ASSERT(complement('T') == 'A', "T -> A");
    ASSERT(complement('C') == 'G', "C -> G");
    ASSERT(complement('G') == 'C', "G -> C");
    ASSERT(complement('N') == 'N', "N -> N");
    ASSERT(complement('X') == 'N', "X -> N (unknown)");
}

// ---- complement_sequence ----
void test_complement_sequence(void) {
    TEST("complement_sequence");

    char dest[64];

    // AACCGGTTAACC -> reverse complement = GGTTAACCGGTT
    complement_sequence("AACCGGTTAACC", dest, 12);
    ASSERT(strcmp(dest, "GGTTAACCGGTT") == 0, "AACCGGTTAACC -> GGTTAACCGGTT");

    // TTGGCCAATTGG -> reverse complement = CCAATTGGCCAA
    complement_sequence("TTGGCCAATTGG", dest, 12);
    ASSERT(strcmp(dest, "CCAATTGGCCAA") == 0, "TTGGCCAATTGG -> CCAATTGGCCAA");

    // Single base
    complement_sequence("A", dest, 1);
    ASSERT(strcmp(dest, "T") == 0, "A -> T (single base)");

    // ATCG -> reverse complement = CGAT
    complement_sequence("ATCG", dest, 4);
    ASSERT(strcmp(dest, "CGAT") == 0, "ATCG -> CGAT");
}

// ---- levenshtein_distance ----
void test_levenshtein_distance(void) {
    TEST("levenshtein_distance");

    // Exact match at start, k=0
    const char *haystack = "AACCGGTTAACCNNNNNN";
    int result = levenshtein_distance(haystack, 18, "AACCGGTTAACC", 12, 0);
    ASSERT(result == 12, "exact match at start returns end position 12");

    // Exact match with offset (10 N's then barcode)
    const char *haystack2 = "NNNNNNNNNNAACCGGTTAACCNNNNN";
    result = levenshtein_distance(haystack2, 26, "AACCGGTTAACC", 12, 0);
    ASSERT(result == 22, "exact match at offset 10 returns 22");

    // 1 mismatch rejected at k=0
    const char *haystack3 = "AACCGTTTAACCNNNNNN";
    result = levenshtein_distance(haystack3, 18, "AACCGGTTAACC", 12, 0);
    ASSERT(result == -1, "1 mismatch rejected at k=0");

    // 1 mismatch accepted at k=1
    result = levenshtein_distance(haystack3, 18, "AACCGGTTAACC", 12, 1);
    ASSERT(result == 12, "1 mismatch accepted at k=1");

    // 2 mismatches at k=2
    const char *haystack4 = "AACCTTTTAACCNNNNNN";
    result = levenshtein_distance(haystack4, 18, "AACCGGTTAACC", 12, 2);
    ASSERT(result == 12, "2 mismatches accepted at k=2");

    // 2 mismatches rejected at k=1
    result = levenshtein_distance(haystack4, 18, "AACCGGTTAACC", 12, 1);
    ASSERT(result == -1, "2 mismatches rejected at k=1");

    // No match at all
    result = levenshtein_distance("NNNNNNNNNNNN", 12, "AACCGGTTAACC", 12, 0);
    ASSERT(result == -1, "no match in all-N haystack at k=0");

    // k > needle_len returns -1
    result = levenshtein_distance("AACCGGTTAACC", 12, "AACCGGTTAACC", 12, 13);
    ASSERT(result == -1, "k > needle_len returns -1");
}

// ---- parse_csv_headers ----
void test_parse_csv_headers(void) {
    TEST("parse_csv_headers");

    int result = parse_csv_headers("./test_barcodes_single.csv");
    ASSERT(result == 1, "single barcode CSV returns 1");

    result = parse_csv_headers("./test_barcodes_dual.csv");
    ASSERT(result == 2, "dual barcode CSV returns 2");

    // Non-existent file
    result = parse_csv_headers("./nonexistent.csv");
    ASSERT(result == -1, "non-existent file returns -1");
}

// ---- is_fastq ----
void test_is_fastq(void) {
    TEST("is_fastq");
    ASSERT(is_fastq("reads.fq") == true, ".fq -> true");
    ASSERT(is_fastq("reads.fastq") == true, ".fastq -> true");
    ASSERT(is_fastq("reads.fasta") == false, ".fasta -> false");
    ASSERT(is_fastq("reads.fq.gz") == true, ".fq.gz -> true");
    ASSERT(is_fastq("reads.fastq.gz") == true, ".fastq.gz -> true");
    ASSERT(is_fastq("reads.txt") == false, ".txt -> false");
}

// ---- average_qual ----
void test_average_qual(void) {
    TEST("average_qual");

    // Uniform quality 'I' = ASCII 73, Phred = 73-33 = 40
    // All same quality -> average should be ~40
    const char *quals = "IIIIIIIIII";
    double q = average_qual(quals, 10);
    ASSERT(fabs(q - 40.0) < 0.01, "uniform quality I -> Phred ~40");

    // Uniform quality '!' = ASCII 33, Phred = 0
    const char *quals2 = "!!!!!!!!!!";
    q = average_qual(quals2, 10);
    ASSERT(fabs(q - 0.0) < 0.01, "uniform quality ! -> Phred ~0");

    // Uniform quality '5' = ASCII 53, Phred = 20
    const char *quals3 = "5555555555";
    q = average_qual(quals3, 10);
    ASSERT(fabs(q - 20.0) < 0.01, "uniform quality 5 -> Phred ~20");
}

// ---- slice ----
void test_slice(void) {
    TEST("slice");

    char dest[64];
    slice("ABCDEFGHIJ", dest, 0, 5);
    ASSERT(strcmp(dest, "ABCDE") == 0, "slice(0,5) -> ABCDE");

    slice("ABCDEFGHIJ", dest, 3, 7);
    ASSERT(strcmp(dest, "DEFG") == 0, "slice(3,7) -> DEFG");

    slice("ABCDEFGHIJ", dest, 0, 10);
    ASSERT(strcmp(dest, "ABCDEFGHIJ") == 0, "slice(0,10) -> full string");

    slice("ABCDEFGHIJ", dest, 9, 10);
    ASSERT(strcmp(dest, "J") == 0, "slice(9,10) -> J");
}

// ---- min ----
void test_min(void) {
    TEST("min");
    ASSERT(min(1, 2, 3) == 1, "min(1,2,3) = 1");
    ASSERT(min(3, 1, 2) == 1, "min(3,1,2) = 1");
    ASSERT(min(2, 3, 1) == 1, "min(2,3,1) = 1");
    ASSERT(min(5, 5, 5) == 5, "min(5,5,5) = 5");
    ASSERT(min(-1, 0, 1) == -1, "min(-1,0,1) = -1");
    ASSERT(min(0, -1, 1) == -1, "min(0,-1,1) = -1");
}

int main(void) {
    printf("=== nanoSweet Unit Tests ===\n\n");

    test_complement();
    test_complement_sequence();
    test_levenshtein_distance();
    test_parse_csv_headers();
    test_is_fastq();
    test_average_qual();
    test_slice();
    test_min();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
