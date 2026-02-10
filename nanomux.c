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

#define FLAG_IMPLEMENTATION
#include "./flag.h"
#include "kseq.h"
#include "thpool.h"
#define COMMON_IMPLEMENTATION
#include "common.h"

#include <zlib.h>
#include <limits.h> 
#include <stdint.h>
#include <pthread.h>

#define READ_BUFFER 10 * 1000
KSEQ_INIT(gzFile, gzread)

typedef struct {
    Barcodes *barcodes;
    Reads *reads;
    size_t start;
    size_t end;
    size_t barcode_pos;
    size_t k;
    bool trim;
    int barcode_schema;
    pthread_mutex_t *bc_mutexes;
} Thread_Data;

void process_reads(void *arg)
{
    Thread_Data *td = (Thread_Data *)arg;
    size_t k = td->k;
    size_t barcode_pos = td->barcode_pos;
    bool trim = td->trim;
    int barcode_schema = td->barcode_schema;

    for (size_t r = td->start; r < td->end; r++) {
        Read *read = &td->reads->items[r];
        const char *first_slice = read->first_slice;
        const char *last_slice = read->last_slice;

        for (size_t bi = 0; bi < td->barcodes->count; bi++) {
            Barcode *b = &td->barcodes->items[bi];
            bool matched = false;

            if (barcode_schema == 1) {
                // Check for barcode in 5' end
                int m5 = levenshtein_distance(first_slice, barcode_pos, b->fw, b->fw_length, k);
                if (m5 != -1) {
                    pthread_mutex_lock(&td->bc_mutexes[bi]);
                    b->counter++;
                    if (trim) {
                        if (!append_read_to_gzip_fastq(b->out_gz, read, m5, read->len)) exit(1);
                    } else {
                        if (!append_read_to_gzip_fastq(b->out_gz, read, 0, read->len)) exit(1);
                    }
                    pthread_mutex_unlock(&td->bc_mutexes[bi]);
                    matched = true;
                } else {
                    // Check for barcode in 3' end
                    int m3 = levenshtein_distance(last_slice, barcode_pos, b->fw_comp, b->fw_length, k);
                    if (m3 != -1) {
                        int slice_end = read->len - barcode_pos + m3 - b->fw_length;
                        if (slice_end > 0) {
                            pthread_mutex_lock(&td->bc_mutexes[bi]);
                            b->counter++;
                            if (trim) {
                                if (!append_read_to_gzip_fastq(b->out_gz, read, 0, slice_end)) exit(1);
                            } else {
                                if (!append_read_to_gzip_fastq(b->out_gz, read, 0, read->len)) exit(1);
                            }
                            pthread_mutex_unlock(&td->bc_mutexes[bi]);
                            matched = true;
                        }
                    }
                }
            } else if (barcode_schema == 2) {
                // fw ------ revcomp(rv)
                int m5fw = levenshtein_distance(first_slice, barcode_pos, b->fw, b->fw_length, k);
                if (m5fw != -1) {
                    int m3rv = levenshtein_distance(last_slice, barcode_pos, b->rv_comp, b->rv_length, k);
                    if (m3rv != -1) {
                        int slice_end = read->len - barcode_pos + m3rv - b->rv_length;
                        if (slice_end > 0) {
                            pthread_mutex_lock(&td->bc_mutexes[bi]);
                            b->counter++;
                            if (trim) {
                                if (!append_read_to_gzip_fastq(b->out_gz, read, m5fw, slice_end)) exit(1);
                            } else {
                                if (!append_read_to_gzip_fastq(b->out_gz, read, 0, read->len)) exit(1);
                            }
                            pthread_mutex_unlock(&td->bc_mutexes[bi]);
                            matched = true;
                        }
                    }
                }
                if (!matched) {
                    // rv ------ revcomp(fw)
                    int m5rv = levenshtein_distance(first_slice, barcode_pos, b->rv, b->rv_length, k);
                    if (m5rv != -1) {
                        int m3fw = levenshtein_distance(last_slice, barcode_pos, b->fw_comp, b->fw_length, k);
                        if (m3fw != -1) {
                            int slice_end = read->len - barcode_pos + m3fw - b->fw_length;
                            if (slice_end > 0) {
                                pthread_mutex_lock(&td->bc_mutexes[bi]);
                                b->counter++;
                                if (trim) {
                                    if (!append_read_to_gzip_fastq(b->out_gz, read, m5rv, slice_end)) exit(1);
                                } else {
                                    if (!append_read_to_gzip_fastq(b->out_gz, read, 0, read->len)) exit(1);
                                }
                                pthread_mutex_unlock(&td->bc_mutexes[bi]);
                                matched = true;
                            }
                        }
                    }
                }
            }
            if (matched) break;
        }
    }
    free(td);
}

void dispatch_reads(threadpool thpool, Barcodes *barcodes, Reads *reads,
                    size_t num_threads, size_t barcode_pos, size_t k,
                    bool trim, int barcode_schema, pthread_mutex_t *bc_mutexes)
{
    if (reads->count == 0) return;
    size_t chunk_size = (reads->count + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; t++) {
        size_t start = t * chunk_size;
        if (start >= reads->count) break;
        size_t end = start + chunk_size;
        if (end > reads->count) end = reads->count;

        Thread_Data *td = malloc(sizeof(Thread_Data));
        if (!td) {
            nob_log(NOB_ERROR, "Failed to allocate thread data");
            exit(1);
        }
        td->barcodes = barcodes;
        td->reads = reads;
        td->start = start;
        td->end = end;
        td->barcode_pos = barcode_pos;
        td->k = k;
        td->trim = trim;
        td->barcode_schema = barcode_schema;
        td->bc_mutexes = bc_mutexes;
        thpool_add_work(thpool, process_reads, (void *)td);
    }
    thpool_wait(thpool);
}

int main(int argc, char **argv) {    

    // flag.h arguments
    char **barcode_file = flag_str("b", "", "Path to barcode file (MANDATORY)");
    char **fastq_file = flag_str("f", "", "Path to fastq file (MANDATORY)");
    char **out_folder = flag_str("o", "", "Name of output folder (MANDATORY)");
    size_t *barcode_pos = flag_size("p", 50, "Position of barcode");
    size_t *k = flag_size("k", 0, "Number of mismatches allowed");
    bool *trim = flag_bool("t", false, "Trim reads from adapters or not");
    size_t *num_threads = flag_size("j", 1, "Number of threads to use");
    bool *help = flag_bool("help", false, "Print this help to stdout and exit with 0");
    bool *version = flag_bool("v", false, "Print the current version");

    if (!flag_parse(argc, argv)) {
        flag_print_options(stderr);
        flag_print_error(stderr);
        return 1;
    }

    if (*help) {
        flag_print_options(stderr);
        return 0;
    }

    if (*version) {
        print_version();
        return 0;
    }

    if (
        strcmp(*barcode_file, "") == 0 || 
        strcmp(*fastq_file, "") == 0 || 
        strcmp(*out_folder, "") == 0
    ) {
        nob_log(NOB_ERROR, "At least one of the mandatory arguments are missing");
        flag_print_options(stderr);
        return 1;
    }

	if (*k >= 4) {
		nob_log(NOB_ERROR, "k cannot be larger than 3");
		return 1;
	}

    nob_log(NOB_INFO, "Running nanomux");
    nob_log(NOB_INFO, "Barcode position: 0 -> %zu", *barcode_pos);
    nob_log(NOB_INFO, "k: %zu", *k);
    const char *trim_option_string = *trim ? "true" : "false";
    nob_log(NOB_INFO, "Trim option: %s", trim_option_string);
    nob_log(NOB_INFO, "threads: %zu", *num_threads);
    printf("\n");

    if (!nob_mkdir_if_not_exists(*out_folder)) {
        nob_log(NOB_ERROR, "exiting");
        return 1;
    }

    
    nob_log(NOB_INFO, "Parsing barcode file %s", *barcode_file);
    
    // ----------------- BARCODES ---------------------------
    int barcode_schema = parse_csv_headers(*barcode_file);
    if (barcode_schema == -1) return 1;
    printf("barcode schema: %d\n", barcode_schema);
    Nob_String_Builder sb = {0};
    Barcodes barcodes = {0};
    if (!parse_barcodes(*barcode_file, &barcodes, &sb, *out_folder)) return 1;
    // validate barcodes
    for (size_t i = 0; i < barcodes.count; i++) {
        Barcode current_bc = barcodes.items[i];
        if (barcode_schema == 2) {
            if (current_bc.rv == NULL) {
                printf("ERROR: Wrong barcode at row: %zu\n", i);
                return 1;
            }
        } else if (barcode_schema == 1)
            if (current_bc.fw == NULL) {
                printf("ERROR: Wrong barcode at row: %zu\n", i);
                return 1;
            }
    }
    
    // ----------------- THREADS ---------------------------
    threadpool thpool = thpool_init(*num_threads);
    if (thpool == NULL) {
        printf("ERROR: Could not init threads\n");
        return 1;
    }

    pthread_mutex_t *bc_mutexes = malloc(barcodes.count * sizeof(pthread_mutex_t));
    if (!bc_mutexes) {
        nob_log(NOB_ERROR, "Failed to allocate mutexes");
        return 1;
    }
    for (size_t i = 0; i < barcodes.count; i++) {
        pthread_mutex_init(&bc_mutexes[i], NULL);
    }
    
    // ----------------- GO THROUGH READS ---------------------------
    gzFile fp = gzopen(*fastq_file, "r"); 
    if (!fp) return 1;
    kseq_t *seq = kseq_init(fp);
    Reads reads = {0};
    int l;
    size_t counter = 0;
    size_t reads_shorter_than_p = 0;

#define REPORT_INTERVAL (1000 * 10)

    while ((l = kseq_read(seq)) >= 0) { 
        counter++;
        if (counter % REPORT_INTERVAL == 0) {
            fprintf(stderr, "\rProcessed: %zu reads", counter);
            fflush(stderr);
        }
        if (seq->seq.l <= *barcode_pos) {
            reads_shorter_than_p++;
            continue;
        }
        Read read = {0};
        read.seq = strdup(seq->seq.s);
        read.name = strdup(seq->name.s);
        read.qual = strdup(seq->qual.s);
        read.len = seq->seq.l;
        read.first_slice = strndup(seq->seq.s, *barcode_pos);
        read.last_slice = strdup(seq->seq.s + seq->seq.l - *barcode_pos);
        
        nob_da_append(&reads, read);
        
        // ----------------- TRIGGER THREADS AND PROCESSING ---------------------------
        if (reads.count >= READ_BUFFER) {
            dispatch_reads(thpool, &barcodes, &reads, *num_threads,
                           *barcode_pos, *k, *trim, barcode_schema, bc_mutexes);

            // Clean up reads
            for (size_t i = 0; i < reads.count; i++) free_read(reads.items[i]);
            reads.count = 0;
        }
    }

    // PROCESS LEFT OVER READS IN BUFFER
    if (reads.count > 0) {
        dispatch_reads(thpool, &barcodes, &reads, *num_threads,
                       *barcode_pos, *k, *trim, barcode_schema, bc_mutexes);
    }
    
    
    // ----------------- LOG TO STDOUT, SUMMARY, MATCHES AND REMOVE EMPTY FILES---------------------------
    FILE *LOG_FILE = open_summary_file(*out_folder, "nanomux.log");
    FILE *S_FILE = open_summary_file(*out_folder, "nanomux_matches.csv");
    
    fprintf(LOG_FILE, "Nanomux\n\n");
    fprintf(LOG_FILE, "Barcodes: %s\n", *barcode_file);
    fprintf(LOG_FILE, "Fastq: %s\n", *fastq_file);
    fprintf(LOG_FILE, "Barcode position: %zu\n", *barcode_pos);
    fprintf(LOG_FILE, "k: %i\n", (int) *k);
    fprintf(LOG_FILE, "Output folder: %s\n", *out_folder);
    fprintf(LOG_FILE, "Trim option: %i\n", *trim);
    printf("\nINFO: Processed %zu reads\n", counter);
    printf("INFO: Reads shorter than p: %zu reads\n", reads_shorter_than_p);
    fprintf(LOG_FILE, "Processed %zu reads\n", counter);
    fprintf(LOG_FILE, "Reads shorter than p: %zu reads\n", reads_shorter_than_p);
    
    fprintf(S_FILE, "barcode,matches\n");
    for (size_t i = 0; i < barcodes.count; i++) {
        size_t bc_count = barcodes.items[i].counter;
        char *bc_name = barcodes.items[i].name;
        fprintf(S_FILE, "%s,%zu\n", bc_name, bc_count);
        printf("%s: %zu\n", bc_name, bc_count);
        // remove file if empty
        if (bc_count == 0) {
            char *bc_file = barcodes.items[i].out_name;
            nob_delete_file(bc_file);
        }
    }
    
    // ----------------- CLEAN-UP ---------------------------
    thpool_destroy(thpool);
    for (size_t i = 0; i < barcodes.count; i++) {
        pthread_mutex_destroy(&bc_mutexes[i]);
    }
    free(bc_mutexes);
    for (size_t i = 0; i < reads.count; i++) free_read(reads.items[i]);
    for (size_t i = 0; i < barcodes.count; i++) {
        gzclose(barcodes.items[i].out_gz);
        free_barcode(&barcodes.items[i]);
    }
    nob_da_free(barcodes);
    nob_da_free(reads);
    fclose(S_FILE);
    fclose(LOG_FILE);
    gzclose(fp);
    
    printf("\n");
    nob_log(NOB_INFO, "nanomux done!\n");
    
    return 0;
}
