// Copyright 2024 William Rosenbaum <william.rosenbaum88@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#define NOB_IMPLEMENTATION
#include "nob.h"
#define FLAG_IMPLEMENTATION
#include "./flag.h"
#include "kseq.h"
#include "thpool.h"

#include <zlib.h>
#include <limits.h> 
#include <stdint.h>
#include <pthread.h>

int min(int a, int b, int c) {
    int min = a;
    if (b < min) {
        min = b;
    }
    if (c < min) {
        min = c;
    }
    return min;
}

// https://stackoverflow.com/questions/8139958/algorithm-to-find-edit-distance-to-all-substrings
int levenshtein_distance(const char *haystack, const char *needle, int k) {
    int needle_len = strlen(needle);
    int haystack_len = strlen(haystack);

    if (k < 0 || k > needle_len) {
        return -1;  
    }

    int dp[needle_len + 1][haystack_len + 1];

    for (int j = 0; j <= haystack_len; j++) {
        dp[0][j] = 0;
    }

    for (int i = 1; i <= needle_len; i++) {
        dp[i][0] = i;
        for (int j = 1; j <= haystack_len; j++) {
            if (needle[i - 1] == haystack[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]);
            }
        }
    }

    for (int j = needle_len; j <= haystack_len; j++) {
        if (dp[needle_len][j] <= k) {
            return j;
        }
    }
    return -1;
}

typedef enum {
    BARCODE_NAME = 0,
    BARCODE_FW,
    BARCODE_RV
} Barcode_Attr;

typedef struct {
    const char *name;
    const char *fw;
    size_t fw_length;
    const char *rv;
    size_t rv_length;
    const char *rv_comp;  
    const char *fw_comp;  
} Barcode;

typedef struct {
    Barcode *items;
    size_t count;
    size_t capacity;
} Barcodes;

typedef struct {
    const char *seq;
    const char *name;
    const char *qual;
} Read;

typedef struct {
    Read *items;
    size_t count;
    size_t capacity;
} Reads;

void free_read(Read *read) {
    if (read) {
        free((void*)read->seq); 
        free((void*)read->qual); 
        free((void*)read->name); 
    }
}

typedef struct {
    Barcode barcode;
    Reads reads;
    const char *output;
    size_t barcode_pos;
    int k;
    int trim;
    FILE *S_FILE;
    pthread_mutex_t *s_mutex;
} NanomuxData;

typedef struct {
    NanomuxData *items;
    size_t count;
    size_t capacity;
} NanomuxDatas;

KSEQ_INIT(gzFile, gzread)
Reads parse_fastq(const char *fastq_file_path) {
    gzFile fp = gzopen(fastq_file_path, "r"); 
    if (!fp) {
        nob_log(NOB_INFO, "Failed to open fastq file, exiting");
        exit(1);
    }
    
    kseq_t *seq = kseq_init(fp); 
    int l;
    Reads reads = {0};

    int counter = 1;
    int print_read_number = 500 * 1000;
    
    while ((l = kseq_read(seq)) >= 0) { 
        Read read = {0};
        read.seq = strdup(seq->seq.s);
        read.name = strdup(seq->name.s);
        read.qual = strdup(seq->qual.s);
        nob_da_append(&reads, read);

        if (counter % print_read_number == 0) nob_log(NOB_INFO, "Parsed: %i reads", counter);
        counter++;
    }
    
    kseq_destroy(seq); 
    gzclose(fp); 
    return reads;
}

char complement(char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N'; 
    }
}

void complement_sequence(const char *src, char *dest, size_t length) {
    for (size_t i = 0; i < length; i++) {
        dest[length - 1 - i] = complement(src[i]);
    }
    dest[length] = '\0';
}

void compute_reverse_complement_rv(Barcode *barcode) {
    char *comp_seq = malloc((barcode->rv_length + 1) * sizeof(char));
    if (comp_seq == NULL) {
        nob_log(NOB_ERROR, "Memory allocation failed for reverse complement, exiting");
        exit(1); 
    }

    complement_sequence(barcode->rv, comp_seq, barcode->rv_length);
    barcode->rv_comp = comp_seq;
}

void compute_reverse_complement_fw(Barcode *barcode) {
    char *comp_seq = malloc((barcode->fw_length + 1) * sizeof(char));
    if (comp_seq == NULL) {
        nob_log(NOB_ERROR, "Memory allocation failed for reverse complement, exiting");
        exit(1); 
    }

    complement_sequence(barcode->fw, comp_seq, barcode->fw_length);
    barcode->fw_comp = comp_seq;
}

Barcodes parse_barcodes(Nob_String_View content) {
    Barcodes barcodes = {0};

    while (content.count > 0) {
        Nob_String_View line = nob_sv_chop_by_delim(&content, '\n');
        Barcode barcode = {0};

        for (int i = 0; line.count > 0; ++i) {
            Nob_String_View attr = nob_sv_chop_by_delim(&line, ',');
            switch (i) {
                case BARCODE_NAME:
                    barcode.name = nob_temp_sv_to_cstr(attr);
                    break;
                case BARCODE_FW:
                    barcode.fw = nob_temp_sv_to_cstr(attr);
                    barcode.fw_length = attr.count;
                    compute_reverse_complement_fw(&barcode); 
                    break;
                case BARCODE_RV:
                    barcode.rv = nob_temp_sv_to_cstr(attr);
                    barcode.rv_length = attr.count;
                    compute_reverse_complement_rv(&barcode); 
                    break;
                default:
                    break;
            }
        }
        nob_da_append(&barcodes, barcode);
    }

    return barcodes;
}

int num_barcode_fields(const char *csv) {
    Nob_String_Builder sb = {0};
    if (!nob_read_entire_file(csv, &sb)) return 1;
    Nob_String_View content = nob_sv_from_parts(sb.items, sb.count);
    Nob_String_View line = nob_sv_chop_by_delim(&content, '\n');
    
    int number_fields = 0;
    while (line.count > 0) {
        nob_sv_chop_by_delim(&line, ',');
        number_fields++;
    }
    
    return number_fields;
}

void slice(const char* src, char* dest, size_t start, size_t end) {
    size_t length = end - start;
    memcpy(dest, src + start, length);
    dest[length] = '\0';  
}

int append_read_to_gzip_fastq(gzFile gzfp, Read *read, int start, int end) {
    int result = 0;
    int length = strlen(read->seq);

    if (start < 0) start = 0;
    if (end > length) end = length;
    size_t trimmed_length = end - start;

    // make sure that we have enough space for the name if the read is short
    size_t buffer_size = trimmed_length + 100; 
    char *buffer = (char *)malloc(buffer_size);

    if (!buffer) {
        nob_log(NOB_ERROR, "Memory allocation failed, buy more RAM lol");
        return 1;
    }

    int len = snprintf(buffer, buffer_size, "@%s\n", read->name);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write name to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%.*s\n", (int)trimmed_length, read->seq + start);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write sequence to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "+\n");
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write separator to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%.*s\n", (int)trimmed_length, read->qual + start);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write quality to gzip file");
        nob_return_defer(1);
    }

defer:
    free(buffer);
    return result;
}

void process_dual_barcode(
    const Barcode b,
    const Reads reads,
    const char *output,
    const size_t barcode_pos,
    const int k,
    const bool trim,
    FILE *s_file,
    pthread_mutex_t *s_mutex
) {
    int counter = 0;
    
    // save fastq
    char fastq_name[512]; 
    snprintf(fastq_name, sizeof(fastq_name), "%s/%s.fq.gz", output, b.name);
    
    gzFile new_fastq = gzopen(fastq_name, "ab");
    if (!new_fastq) {
        nob_log(NOB_INFO, "failed to open new fastq file for appending, exiting");
        exit(1);
    }
    
    char target_slice[barcode_pos + 10];
    
    for (size_t j = 0; j < reads.count; ++j) {
        Read r = reads.items[j];
        int length = strlen(r.seq);
        
        // fw ------ revcomp(rv)
        // fw
        slice(r.seq, target_slice, 0, barcode_pos);
        int match_first_fw = levenshtein_distance(target_slice, b.fw, k); 
        if (match_first_fw != -1) {
            slice(r.seq, target_slice, length - barcode_pos, length);
            // revcomp(rv)
            int match_last_fw = levenshtein_distance(target_slice, b.rv_comp, k);
            if (match_last_fw != -1) {
                counter++;
                int slice_end = length - barcode_pos + match_last_fw - b.rv_length;
                if (trim) {
                    append_read_to_gzip_fastq(new_fastq, &r, match_first_fw, slice_end);
                } else {
                    append_read_to_gzip_fastq(new_fastq, &r, 0, length);
                }
                continue;
            }
        }
        
        // rv ------ revcomp(fw)
        slice(r.seq, target_slice, 0, barcode_pos);
        int match_first_rv = levenshtein_distance(target_slice, b.rv, k); 
        if (match_first_rv != -1) {
            slice(r.seq, target_slice, length - barcode_pos, length);
            // revcomp(rv)
            int match_last_rv = levenshtein_distance(target_slice, b.fw_comp, k);
            if (match_last_rv != -1) {
                counter++;
                int slice_end = length - barcode_pos + match_last_rv - b.fw_length;
                if (trim) {
                    append_read_to_gzip_fastq(new_fastq, &r, match_first_rv, slice_end);
                } else {
                    append_read_to_gzip_fastq(new_fastq, &r, 0, length);
                }
            }
        }
    }
        
    gzclose(new_fastq);
    
    pthread_mutex_lock(s_mutex);
    nob_log(NOB_INFO, "barcode: %30s matches: %i", b.name, counter);
    fprintf(s_file, "%s,%i\n", b.name, counter);
    pthread_mutex_unlock(s_mutex);
}

void process_single_barcode(
    const Barcode b,
    const Reads reads,
    const char *output,
    const size_t barcode_pos,
    const int k,
    const bool trim,
    FILE *s_file,
    pthread_mutex_t *s_mutex
) {
    int counter = 0;
    
    // save fastq
    char fastq_name[512]; 
    snprintf(fastq_name, sizeof(fastq_name), "%s/%s.fq.gz", output, b.name);
    
    gzFile new_fastq = gzopen(fastq_name, "ab");
    if (!new_fastq) {
        nob_log(NOB_INFO, "failed to open new fastq_file for appending, exiting");
        exit(1);
    }
    
    char target_slice[barcode_pos + 10];
    
    for (size_t j = 0; j < reads.count; ++j) {
        Read r = reads.items[j];
        int length = strlen(r.seq);
        
        // Check for barcode in 5' end
        slice(r.seq, target_slice, 0, barcode_pos);
        int match_first_fw = levenshtein_distance(target_slice, b.fw, k);
        if (match_first_fw != -1) {
            counter++;
            if (trim) {
                append_read_to_gzip_fastq(new_fastq, &r, match_first_fw, length);
            } else {
                append_read_to_gzip_fastq(new_fastq, &r, 0, length);
            }
            continue;
        }
        
        // Check for barcode in 3' end
        slice(r.seq, target_slice, length - barcode_pos, length);
        int match_last_rv = levenshtein_distance(target_slice, b.fw_comp, k);
        if (match_last_rv != -1) {
            counter++;
            int slice_end = length - barcode_pos + match_last_rv - b.fw_length;
            if (trim) {
                append_read_to_gzip_fastq(new_fastq, &r, 0, slice_end);
            } else {
                append_read_to_gzip_fastq(new_fastq, &r, 0, length);
            }
        }
    }
    
    gzclose(new_fastq);
    
    pthread_mutex_lock(s_mutex);
    nob_log(NOB_INFO, "barcode: %30s matches: %i", b.name, counter);
    fprintf(s_file, "%s,%i\n", b.name, counter);
    pthread_mutex_unlock(s_mutex);
}

void run_nanomux_dual(void *arg) {
    NanomuxData *data = (NanomuxData *)arg;
    Reads local_reads = data->reads;

    process_dual_barcode(
        data->barcode, 
        local_reads, 
        data->output, 
        data->barcode_pos,
        data->k, 
        data->trim,
        data->S_FILE,
        data->s_mutex
    );
}

void run_nanomux_single(void *arg) {
    NanomuxData *data = (NanomuxData *)arg;
    Reads local_reads = data->reads;

    process_single_barcode(
        data->barcode, 
        local_reads, 
        data->output, 
        data->barcode_pos,
        data->k, 
        data->trim,
        data->S_FILE,
        data->s_mutex
    );
}

void print_barcodes_documentation(void) {
    printf("\nYOUR BARCODE FILE IS WRONG!\n");
    printf("You can use either single barcodes, or dual barcodes\n\n");
    printf("Dual barcodes example:\n");
    printf("barcode1,ACTATCTACTA,GAGCATGTCGTA\n");
    printf("barcode2,AGCGTATGCTGGTA,AGCATGCTATCG\n\n");
    printf("Single barcode example:\n");
    printf("barcode1,ACTATCTACTA\n");
    printf("barcode2,AGCGTATGCTGGTA\n");
}

#define LOG_FILE_CAP 512

int main(int argc, char **argv) {    

    // flag.h arguments
    char **barcode_file = flag_str("b", "", "Path to barcode file (MANDATORY)");
    char **fastq_file = flag_str("f", "", "Path to fastq file (MANDATORY)");
    char **output = flag_str("o", "", "Name of output folder (MANDATORY)");
    size_t *barcode_pos = flag_size("p", 50, "Position of barcode");
    size_t *k = flag_size("k", 0, "Number of mismatches allowed");
    bool *trim = flag_bool("t", true, "Trim reads from adapters or not");
    size_t *num_threads = flag_size("j", 1, "Number of threads to use");
    bool *help = flag_bool("help", false, "Print this help to stdout and exit with 0");

    if (!flag_parse(argc, argv)) {
        flag_print_options(stderr);
        flag_print_error(stderr);
        return 1;
    }

    if (*help) {
        flag_print_options(stderr);
        return 0;
    }

    if (strcmp(*barcode_file, "") == 0 || strcmp(*fastq_file, "") == 0 || strcmp(*output, "") == 0) {
        nob_log(NOB_ERROR, "At least one of the mandatory arguments are missing");
        flag_print_options(stderr);
        return 1;
    }

    nob_log(NOB_INFO, "Running nanomux");
    nob_log(NOB_INFO, "Barcode position: 0 -> %zu", *barcode_pos);
    nob_log(NOB_INFO, "k: %i", *k);
    nob_log(NOB_INFO, "Trim option: %i", *trim);
    nob_log(NOB_INFO, "threads: %i", *num_threads);
    printf("\n");

    if (!nob_mkdir_if_not_exists(*output)) {
        nob_log(NOB_ERROR, "exiting");
        return 1;
    }

    char log_file[LOG_FILE_CAP];
    snprintf(log_file, sizeof(log_file), "%s/nanomux.log", *output);
    FILE *LOG_FILE = fopen(log_file, "ab");
    if (LOG_FILE == NULL) {
        nob_log(NOB_ERROR, "Could NOT create log file");
        return 1;
    }
    
    char summary_file[LOG_FILE_CAP];
    snprintf(summary_file, sizeof(summary_file), "%s/nanomux_matches.csv", *output);
    FILE *S_FILE = fopen(summary_file, "ab");
    if (S_FILE == NULL) {
        nob_log(NOB_ERROR, "Could NOT create summary file");
        fclose(LOG_FILE);
        return 1;
    }
    
    fprintf(S_FILE, "barcode,matches\n");
    fprintf(LOG_FILE, "Nanomux\n\n");
    fprintf(LOG_FILE, "Barcodes: %s\n", *barcode_file);
    fprintf(LOG_FILE, "Fastq: %s\n", *fastq_file);
    fprintf(LOG_FILE, "Barcode position: %zu\n", *barcode_pos);
    fprintf(LOG_FILE, "k: %i\n", (int) *k);
    fprintf(LOG_FILE, "Output folder: %s\n", *output);
    fprintf(LOG_FILE, "Trim option: %i\n", *trim);
    
    nob_log(NOB_INFO, "Parsing barcode file %s", *barcode_file);
    Nob_String_Builder sb = {0};
    if (!nob_read_entire_file(*barcode_file, &sb)) return 1;
    Nob_String_View content = nob_sv_from_parts(sb.items, sb.count);
    Barcodes barcodes = parse_barcodes(content);
    
    // If the number fields in barcode file is 3 --> the file contains fw and rv
    // If the number fields in barcode file is 2 --> the file contains only fw
    char *barcode_type;
    int num_fields = num_barcode_fields(*barcode_file);
    if (num_fields == 2) {
        barcode_type = "single";
        nob_log(NOB_INFO, "Barcode file contains single barcode\n");
    } else if (num_fields == 3) {
        barcode_type = "dual";
        nob_log(NOB_INFO, "Barcode file contains dual barcodes\n");
    } else {
        print_barcodes_documentation();
        return 1;
    }
    
    nob_log(NOB_INFO, "Parsing fastq file %s", *fastq_file);
    Reads reads = parse_fastq(*fastq_file);
    nob_log(NOB_INFO, "Number of reads: %zu", reads.count);
    fprintf(LOG_FILE, "Number of reads: %zu\n", reads.count);
    printf("\n");

    nob_log(NOB_INFO, "Generating threadpool with %i threads", *num_threads);
    printf("\n");
    threadpool thpool = thpool_init(*num_threads);
    
    NanomuxDatas nanomux_datas = {0};
    pthread_mutex_t s_mutex = PTHREAD_MUTEX_INITIALIZER;
    for (size_t i = 0; i < barcodes.count; i++) {
        NanomuxData nanomux_data = {
            .barcode = barcodes.items[i],
            .reads = reads,
            .output = *output,
            .barcode_pos = *barcode_pos,
            .k = (int)*k,
            .trim = *trim,
            .S_FILE = S_FILE,
            .s_mutex = &s_mutex,
        };
        nob_da_append(&nanomux_datas, nanomux_data);
    }
    
    if (strcmp(barcode_type, "dual") == 0) {
        for (size_t i = 0; i < nanomux_datas.count; i++) {
        	thpool_add_work(thpool, run_nanomux_dual, (void *)&nanomux_datas.items[i]);
        }
    } else if (strcmp(barcode_type, "single") == 0) {
        for (size_t i = 0; i < nanomux_datas.count; i++) {
        	thpool_add_work(thpool, run_nanomux_single, (void *)&nanomux_datas.items[i]);
        }
    } else {
        // unreachable
        print_barcodes_documentation();
        return 1;
    }

    thpool_wait(thpool);
    thpool_destroy(thpool);

    pthread_mutex_destroy(&s_mutex);

    nob_da_free(barcodes);
    nob_da_free(reads);
    nob_da_free(nanomux_datas);

    fclose(S_FILE);
    fclose(LOG_FILE);

    printf("\n");
    nob_log(NOB_INFO, "nanomux done!\n");
    
    return 0;
}