#ifndef COMMON_H_
#define COMMON_H_

#define NOB_IMPLEMENTATION
#include "nob.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <pthread.h>
#include <math.h>


// ----------------------------------------------------------------------------
#define FILE_CAP 524
#define VERSION "2.0.0"

typedef struct {
    char *name;

    char *fw;
    char *fw_comp;  
    size_t fw_length;

    char *rv;
    char *rv_comp;  
    size_t rv_length;

    char out_name[512];
    gzFile out_gz;

    size_t counter;
} Barcode;

typedef struct{
    Barcode *items;
    size_t count;
    size_t capacity;
} Barcodes;

typedef struct {
    const char *seq;
    const char *name;
    const char *qual;
    const char *first_slice;
    const char *last_slice;
    size_t len;
} Read;

typedef struct {
    Read *items;
    size_t count;
    size_t capacity;
} Reads;

bool append_read_to_gzip_fastq(gzFile gzfp, Read *read, int start, int end);
void print_barcode_documentation(void);
void slice_str(const char * str, char * buffer, size_t start, size_t end);
void slice(const char* src, char* dest, size_t start, size_t end);
char complement(const char nucleotide);
void complement_sequence(char *src, char *dest, size_t length);
bool parse_barcodes(const char *bc_path, Barcodes *barcodes, Nob_String_Builder *sb, char *outdir);
int parse_csv_headers(const char *barcode_path);
void close_gz_files(Barcode *bc);
void free_barcode(Barcode *bc);
static inline int min(int a, int b, int c);
int levenshtein_distance(const char *haystack, size_t haystack_len, const char *needle, size_t needle_len, size_t k);
FILE* open_summary_file(const char *out_folder, const char *filename);
void free_read(Read read);
char *basename(char const *path);
double average_qual(const char *quals, size_t len);
bool is_fastq(const char *file);
bool must_be_digit(const char *arg);
void print_version(void);

#endif // COMMON_H_

// ----------------------------------------------------------------------------

#ifdef COMMON_IMPLEMENTATION

char complement(const char nucleotide) 
{
    switch (nucleotide) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N'; 
    }
}

void complement_sequence(char *src, char *dest, size_t length) 
{
    for (size_t i = 0; i < length; i++) {
        dest[length - 1 - i] = complement(src[i]);
    }
    dest[length] = '\0';
}

bool parse_barcodes(const char *bc_path, Barcodes *barcodes, Nob_String_Builder *sb, char *outdir)
{
    if (!nob_read_entire_file(bc_path, sb)) return false;

    Nob_String_View content = nob_sb_to_sv(*sb);

    bool first_line = true;  
    while(content.count > 0) {
        Nob_String_View line = nob_sv_chop_by_delim(&content, '\n');
        
        // skip the header
        if (first_line) {
            first_line = false;
            continue;
        }
        
        Barcode barcode = {0};
        
        for (size_t i = 0; i < line.count; i++) {
            Nob_String_View field = nob_sv_chop_by_delim(&line, ',');
            const char *field_cstr = nob_temp_sv_to_cstr(field);
            size_t field_len = strlen(field_cstr);
            
            switch (i) {
                case 0:
                    barcode.name = strdup(field_cstr);
                    break;
                case 1:
                    barcode.fw_length = field_len;
                    barcode.fw = strdup(field_cstr);
                    barcode.fw_comp = malloc(field_len + 1);
                    if (barcode.fw_comp == NULL) {
                        free(barcode.name);
                        free(barcode.fw);
                        return false;
                    }
                    complement_sequence(barcode.fw, barcode.fw_comp, barcode.fw_length);
                    break;
                case 2:
                    barcode.rv_length = strlen(field_cstr);
                    barcode.rv = strdup(field_cstr);
                    barcode.rv_comp = malloc(field_len + 1);
                    if (barcode.rv_comp == NULL) {
                        free(barcode.name);
                        free(barcode.fw);
                        free(barcode.fw_comp);
                        free(barcode.rv);
                        return false;
                    }
                    complement_sequence(barcode.rv, barcode.rv_comp, barcode.rv_length);
                    break;
                default: 
                    printf("ERROR: your barcodes contains %zu fields. It should be 3.\n", i + 1);
                    print_barcode_documentation();
                    return false;
            }
        }
        // add the new gz file to write to later.
        snprintf(barcode.out_name, sizeof(barcode.out_name), "%s/%s.fq.gz", outdir, barcode.name);
        barcode.out_gz = gzopen(barcode.out_name, "ab");
        if (!barcode.out_gz) {
            printf("ERROR: Could not open %s to write to\n", barcode.out_name);
            return false;
        }
        nob_da_append(barcodes, barcode);
    }
    return true;
}

void close_gz_files(Barcode *bc)
{
    if (bc->out_gz) gzclose(bc->out_gz);
}

void free_barcode(Barcode *bc)
{
    free(bc->name);
    free(bc->fw);
    free(bc->fw_comp);
    free(bc->rv);
    free(bc->rv_comp);
}

void print_barcode_documentation(void) 
{
    printf("You can use either single barcodes, or dual barcodes\n\n");
    printf("Dual barcodes example:\n");
    printf("name,forward,reverse\n");
    printf("barcode1,ACTATCTACTA,GAGCATGTCGTA\n");
    printf("barcode2,AGCGTATGCTGGTA,AGCATGCTATCG\n\n");
    printf("Single barcode example:\n");
    printf("name,forward\n");
    printf("barcode1,ACTATCTACTA\n");
    printf("barcode2,AGCGTATGCTGGTA\n");
}

// https://stackoverflow.com/questions/8139958/algorithm-to-find-edit-distance-to-all-substrings
int levenshtein_distance(const char *haystack, size_t haystack_len, const char *needle, size_t needle_len, size_t k) 
{
    if (k > needle_len) return -1;
    
    size_t dp[needle_len + 1][haystack_len + 1];

    for (size_t j = 0; j <= haystack_len; j++) {
        dp[0][j] = 0;
    }

    for (size_t i = 1; i <= needle_len; i++) {
        dp[i][0] = i;
        for (size_t j = 1; j <= haystack_len; j++) {
            if (needle[i - 1] == haystack[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]);
            }
        }
    }

    for (size_t j = needle_len; j <= haystack_len; j++) {
        if (dp[needle_len][j] <= k) return (int)j;
    }
    return -1;
}

static inline int min(int a, int b, int c) 
{
    int min = a;
    if (b < min) {
        min = b;
    }
    if (c < min) {
        min = c;
    }
    return min;
}

int parse_csv_headers(const char *barcode_path) 
{
    FILE *f = fopen(barcode_path, "r");
    if (!f) {
        printf("Error: couldn't open file\n");
        return -1;
    }
    
    char cols[10][100];  
    int count = 0;
    int c;
    
    while (fscanf(f, "%99[^,\n]", cols[count]) == 1) {
        count++;
        if (count >= 10) break;
        c = fgetc(f);
        if (c == '\n' || c == EOF) break;
    }
    
    fclose(f);

    if (!(count == 2 || count == 3)) {
        printf("ERROR: Wrong amount of barcode headers\n");
        print_barcode_documentation();
        return -1;
    }

    // single barcode
    if (count == 2) {
        if (strcmp(cols[0], "name") != 0 || strcmp(cols[1], "forward") != 0) {
            printf("ERROR: Headers of barcode file have incorrect headers\n");
            print_barcode_documentation();
            return -1;
        }
        return 1;
    // dual barcode
    } else {
        if (strcmp(cols[0], "name") != 0 || strcmp(cols[1], "forward") != 0 || strcmp(cols[2], "reverse") != 0) {
            printf("ERROR: Headers of barcode file have incorrect headers\n");
            print_barcode_documentation();
            return -1;
        }
       return 2; 
    }
    //UNREACHABLE
}

void slice(const char* src, char* dest, size_t start, size_t end) 
{
    size_t length = end - start;
    memcpy(dest, src + start, length);
    dest[length] = '\0';  
}

FILE* open_summary_file(const char *out_folder, const char *filename) 
{
    char summary_file[FILE_CAP];
    snprintf(summary_file, sizeof(summary_file), "%s/%s", out_folder, filename);
    
    FILE *S_FILE = fopen(summary_file, "ab");
    if (S_FILE == NULL) {
        nob_log(NOB_ERROR, "Could NOT create summary file");
        return NULL;
    }
    
    return S_FILE;
}


void free_read(Read read) 
{
    free((void *)read.seq);
    free((void *)read.name);
    free((void *)read.qual);
    free((void *)read.first_slice);
    free((void *)read.last_slice);
}

bool append_read_to_gzip_fastq(gzFile gzfp, Read *read, int start, int end) 
{
    int length = read->len;  
    if (start < 0) start = 0;
    if (end > length) end = length;
    if (start >= end) {
        printf("ERROR: Invalid trim range: start=%d, end=%d\n", start, end);
        return false;
    }
    
    size_t trimmed_length = end - start;
    int ret = gzprintf(
        gzfp, 
        "@%s\n"
        "%.*s\n"
        "+\n"
        "%.*s\n", 
        read->name, 
        (int)trimmed_length, read->seq + start,
        (int)trimmed_length, read->qual + start
    );
    
    if (ret < 0) {
        printf("ERROR: Failed to write FASTQ record\n");
        return false;
    }
    return true;
}

char *basename(char const *path) 
{
    char *s = strrchr(path, '/');
    if (!s)
        return strdup(path);
    else
        return strdup(s + 1);
}

double average_qual(const char *quals, size_t len) 
{
    double probability_sum = 0.0;
    for (size_t i = 0; i < len; i++) {
        int phred_score = quals[i] - 33;  
        probability_sum += pow(10.0, phred_score / -10.0);
    }
    return log10(probability_sum / len) * -10.0;
}

bool is_fastq(const char *file) 
{
    return strstr(file, "fastq") || strstr(file, "fq");
}

bool must_be_digit(const char *arg) 
{
    for(; *arg != '\0'; arg++) {
        if (!isdigit(*arg)) return false;
    }
    return true;
}

void print_version(void)
{
    printf("v. %s\n", VERSION);
}


#endif // COMMON_IMPLEMENTATION