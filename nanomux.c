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
#include <sys/resource.h>

#define READ_BUFFER 10 * 1000
#define FILE_DESCRIPTOR_RESERVE 16
KSEQ_INIT(gzFile, gzread)

typedef struct {
    Read *read;
    int start;
    int end;
} Pending_Write;

typedef struct {
    Pending_Write *items;
    size_t count;
    size_t capacity;
} Pending_Writes;

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t available;
    size_t active;
    size_t limit;
} Output_Gate;

typedef struct {
    Barcode *barcode;
    Reads *reads;
    size_t barcode_pos;
    size_t k;
    bool trim;
    int barcode_schema;
    pthread_mutex_t *error_mutex;
    bool *worker_failed;
    Output_Gate *output_gate;
} Thread_Data;

static void mark_worker_failed(Thread_Data *td)
{
    pthread_mutex_lock(td->error_mutex);
    *td->worker_failed = true;
    pthread_mutex_unlock(td->error_mutex);
}

static bool output_gate_init(Output_Gate *gate, size_t limit)
{
    int status = pthread_mutex_init(&gate->mutex, NULL);
    if (status != 0) {
        nob_log(NOB_ERROR, "Could not initialize output mutex: %s", strerror(status));
        return false;
    }

    status = pthread_cond_init(&gate->available, NULL);
    if (status != 0) {
        nob_log(NOB_ERROR, "Could not initialize output condition: %s", strerror(status));
        pthread_mutex_destroy(&gate->mutex);
        return false;
    }

    gate->active = 0;
    gate->limit = limit;
    return true;
}

static bool output_gate_destroy(Output_Gate *gate)
{
    int cond_status = pthread_cond_destroy(&gate->available);
    int mutex_status = pthread_mutex_destroy(&gate->mutex);
    if (cond_status != 0 || mutex_status != 0) {
        nob_log(NOB_ERROR, "Could not destroy output synchronization primitives");
        return false;
    }
    return true;
}

static void output_gate_acquire(Output_Gate *gate)
{
    pthread_mutex_lock(&gate->mutex);
    while (gate->active >= gate->limit) {
        pthread_cond_wait(&gate->available, &gate->mutex);
    }
    gate->active++;
    pthread_mutex_unlock(&gate->mutex);
}

static void output_gate_release(Output_Gate *gate)
{
    pthread_mutex_lock(&gate->mutex);
    gate->active--;
    pthread_cond_signal(&gate->available);
    pthread_mutex_unlock(&gate->mutex);
}

static void queue_barcode_read(Pending_Writes *writes, Read *read, int start, int end)
{
    Pending_Write pending = {
        .read = read,
        .start = start,
        .end = end,
    };
    nob_da_append(writes, pending);
}

static size_t safe_output_handle_count(size_t requested_threads)
{
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0 || limit.rlim_cur == RLIM_INFINITY) {
        return requested_threads;
    }

    if (limit.rlim_cur <= FILE_DESCRIPTOR_RESERVE) return 1;

    rlim_t available = limit.rlim_cur - FILE_DESCRIPTOR_RESERVE;
    if ((rlim_t)requested_threads > available) return (size_t)available;
    return requested_threads;
}

static bool write_pending_reads(Thread_Data *td, Pending_Writes *writes)
{
    if (writes->count == 0) return true;

    bool success = true;
    gzFile output = NULL;
    output_gate_acquire(td->output_gate);

    errno = 0;
    output = gzopen(td->barcode->out_name, "ab");
    if (output == NULL) {
        if (errno != 0) {
            nob_log(NOB_ERROR, "Could not open %s for writing: %s", td->barcode->out_name, strerror(errno));
        } else {
            nob_log(NOB_ERROR, "Could not open %s for writing", td->barcode->out_name);
        }
        success = false;
        goto cleanup;
    }

    for (size_t i = 0; i < writes->count; i++) {
        Pending_Write *pending = &writes->items[i];
        if (!append_read_to_gzip_fastq(output, pending->read, pending->start, pending->end)) {
            success = false;
            break;
        }
    }

cleanup:
    if (output != NULL) {
        int close_status = gzclose(output);
        if (close_status != Z_OK) {
            nob_log(NOB_ERROR, "Failed to close %s: %s", td->barcode->out_name, zError(close_status));
            success = false;
        }
    }
    output_gate_release(td->output_gate);
    return success;
}

void process_barcode(void *arg) 
{
    Thread_Data *td = (Thread_Data *)arg;
    Barcode *b = td->barcode;
    size_t k = td->k;
    size_t barcode_pos = td->barcode_pos;
    bool trim = td->trim;
    int barcode_schema = td->barcode_schema;
    bool success = true;
    Pending_Writes pending_writes = {0};

    // Single barcode processing
    if (!(barcode_schema == 1 || barcode_schema == 2)) {
        nob_log(NOB_ERROR, "Wrong barcode schema");
        success = false;
        goto cleanup;
    } 
    if (barcode_schema == 1) {
        for (size_t i = 0; i < td->reads->count; i++) {
            Read *read = &td->reads->items[i];
            const char *first_read_slice = read->first_slice;
            const char *last_read_slice = read->last_slice;
            
            // Check for barcode in 5' end
            int match_first_fw = levenshtein_distance(first_read_slice, barcode_pos, b->fw, b->fw_length, k);
            if (match_first_fw != -1) {
                b->counter++;
                if (trim) {
                    queue_barcode_read(&pending_writes, read, match_first_fw, read->len);
                } else {
                    queue_barcode_read(&pending_writes, read, 0, read->len);
                }
            } else {
                // Check for barcode in 3' end
                int match_last_rv = levenshtein_distance(last_read_slice, barcode_pos, b->fw_comp, b->fw_length, k);
                if (match_last_rv != -1) {
                    b->counter++;
                    int slice_end = read->len - barcode_pos + match_last_rv - b->fw_length;
                    if (slice_end <= 0) continue;
                    if (trim) {
                        queue_barcode_read(&pending_writes, read, 0, slice_end);
                    } else {
                        queue_barcode_read(&pending_writes, read, 0, read->len);
                    }
                }
            }
        }
    }

    // Dual barcode processing
    if (barcode_schema == 2) {
        for (size_t i = 0; i < td->reads->count; i++) {
            Read *read = &td->reads->items[i];
            const char *first_read_slice = read->first_slice;
            const char *last_read_slice = read->last_slice;
            
            // fw ------ revcomp(rv)
            int match_first_fw = levenshtein_distance(first_read_slice, barcode_pos, b->fw, b->fw_length, k); 
            //printf("DEBUG levenshtein: haystack='%s' (len %zu), needle='%s' (len %zu), k=%zu, result=%d\n", first_read_slice, barcode_pos, b->fw, b->fw_length, k, match_first_fw);
            if (match_first_fw != -1) {
                // revcomp(rv)
                int match_last_fw = levenshtein_distance(last_read_slice, barcode_pos, b->rv_comp, b->rv_length, k);
                if (match_last_fw != -1) {
                    b->counter++;
                    int slice_end = read->len - barcode_pos + match_last_fw - b->rv_length;
                    if (slice_end <= 0) continue;
                    if (trim) {
                        queue_barcode_read(&pending_writes, read, match_first_fw, slice_end);
                    } else {
                        queue_barcode_read(&pending_writes, read, 0, read->len);
                    }
                }
            } else {
                // rv ------ revcomp(fw)
                int match_first_rv = levenshtein_distance(first_read_slice, barcode_pos, b->rv, b->rv_length, k);
                if (match_first_rv != -1) {
                    // revcomp(fw)
                    int match_last_rv = levenshtein_distance(last_read_slice, barcode_pos, b->fw_comp, b->fw_length, k);
                    if (match_last_rv != -1) {
                        b->counter++;
                        int slice_end = read->len - barcode_pos + match_last_rv - b->fw_length;
                        if (slice_end <= 0) continue;
                        if (trim) {
                            queue_barcode_read(&pending_writes, read, match_first_rv, slice_end);
                        } else {
                            queue_barcode_read(&pending_writes, read, 0, read->len);
                        }
                    }
                }
            }
        }
    }

cleanup:
    if (success && !write_pending_reads(td, &pending_writes)) success = false;
    nob_da_free(pending_writes);

    if (!success) {
        mark_worker_failed(td);
    }
    free(td);
}

static bool dispatch_barcodes(threadpool thpool, Barcodes *barcodes, Reads *reads,
                              size_t barcode_pos, size_t k, bool trim,
                              int barcode_schema, pthread_mutex_t *error_mutex,
                              bool *worker_failed, Output_Gate *output_gate)
{
    bool all_jobs_scheduled = true;

    for (size_t i = 0; i < barcodes->count; i++) {
        Thread_Data *td = malloc(sizeof(Thread_Data));
        if (td == NULL) {
            nob_log(NOB_ERROR, "Failed to allocate thread data");
            all_jobs_scheduled = false;
            break;
        }

        td->barcode = &barcodes->items[i];
        td->reads = reads;
        td->barcode_pos = barcode_pos;
        td->k = k;
        td->trim = trim;
        td->barcode_schema = barcode_schema;
        td->error_mutex = error_mutex;
        td->worker_failed = worker_failed;
        td->output_gate = output_gate;

        if (thpool_add_work(thpool, process_barcode, (void *)td) != 0) {
            nob_log(NOB_ERROR, "Failed to schedule barcode processing");
            free(td);
            all_jobs_scheduled = false;
            break;
        }
    }

    thpool_wait(thpool);

    pthread_mutex_lock(error_mutex);
    bool failed = *worker_failed;
    pthread_mutex_unlock(error_mutex);
    return all_jobs_scheduled && !failed;
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

    if (*num_threads == 0) {
        nob_log(NOB_ERROR, "Number of threads must be at least 1");
        return 1;
    }

    int result = 1;
    int barcode_schema = -1;
    Nob_String_Builder sb = {0};
    Barcodes barcodes = {0};
    Reads reads = {0};
    threadpool thpool = NULL;
    gzFile fp = NULL;
    kseq_t *seq = NULL;
    FILE *LOG_FILE = NULL;
    FILE *S_FILE = NULL;
    pthread_mutex_t error_mutex;
    Output_Gate output_gate;
    bool error_mutex_initialized = false;
    bool output_gate_initialized = false;
    bool worker_failed = false;
    size_t counter = 0;
    size_t reads_shorter_than_p = 0;
    size_t output_handle_count = safe_output_handle_count(*num_threads);

    nob_log(NOB_INFO, "Running nanomux");
    nob_log(NOB_INFO, "Barcode position: 0 -> %zu", *barcode_pos);
    nob_log(NOB_INFO, "k: %zu", *k);
    const char *trim_option_string = *trim ? "true" : "false";
    nob_log(NOB_INFO, "Trim option: %s", trim_option_string);
    if (output_handle_count != *num_threads) {
        nob_log(NOB_WARNING, "Concurrent gzip outputs limited to %zu by the file descriptor limit", output_handle_count);
    }
    nob_log(NOB_INFO, "threads: %zu", *num_threads);
    printf("\n");

    if (!nob_mkdir_if_not_exists(*out_folder)) {
        nob_log(NOB_ERROR, "exiting");
        goto cleanup;
    }

    
    nob_log(NOB_INFO, "Parsing barcode file %s", *barcode_file);
    
    // ----------------- BARCODES ---------------------------
    barcode_schema = parse_csv_headers(*barcode_file);
    if (barcode_schema == -1) goto cleanup;
    printf("barcode schema: %d\n", barcode_schema);
    if (!parse_barcodes(*barcode_file, &barcodes, &sb, *out_folder)) goto cleanup;
    // validate barcodes
    for (size_t i = 0; i < barcodes.count; i++) {
        Barcode current_bc = barcodes.items[i];
        if (barcode_schema == 2) {
            if (current_bc.rv == NULL) {
                printf("ERROR: Wrong barcode at row: %zu\n", i);
                goto cleanup;
            }
        } else if (barcode_schema == 1)
            if (current_bc.fw == NULL) {
                printf("ERROR: Wrong barcode at row: %zu\n", i);
                goto cleanup;
            }
    }
    
    // ----------------- THREADS ---------------------------
    thpool = thpool_init(*num_threads);
    if (thpool == NULL) {
        nob_log(NOB_ERROR, "Could not initialize threads");
        goto cleanup;
    }

    int mutex_status = pthread_mutex_init(&error_mutex, NULL);
    if (mutex_status != 0) {
        nob_log(NOB_ERROR, "Could not initialize worker error mutex: %s", strerror(mutex_status));
        goto cleanup;
    }
    error_mutex_initialized = true;

    if (!output_gate_init(&output_gate, output_handle_count)) {
        goto cleanup;
    }
    output_gate_initialized = true;
    
    // ----------------- GO THROUGH READS ---------------------------
    errno = 0;
    fp = gzopen(*fastq_file, "r");
    if (fp == NULL) {
        nob_log(NOB_ERROR, "Could not open %s for reading: %s", *fastq_file, strerror(errno));
        goto cleanup;
    }
    seq = kseq_init(fp);
    if (seq == NULL) {
        nob_log(NOB_ERROR, "Could not initialize FASTQ reader");
        goto cleanup;
    }
    int l;

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
            if (!dispatch_barcodes(thpool, &barcodes, &reads, *barcode_pos, *k,
                                   *trim, barcode_schema, &error_mutex, &worker_failed, &output_gate)) {
                goto cleanup;
            }
            
            // Clean up reads
            for (size_t i = 0; i < reads.count; i++) free_read(reads.items[i]);
            reads.count = 0;
        }
    }

    if (l < -1) {
        nob_log(NOB_ERROR, "Failed while reading FASTQ input");
        goto cleanup;
    }

    // PROCESS LEFT OVER READS IN BUFFER
    if (reads.count > 0) {
        if (!dispatch_barcodes(thpool, &barcodes, &reads, *barcode_pos, *k,
                               *trim, barcode_schema, &error_mutex, &worker_failed, &output_gate)) {
            goto cleanup;
        }
    }
    
    
    // ----------------- LOG TO STDOUT, SUMMARY, MATCHES AND REMOVE EMPTY FILES---------------------------
    LOG_FILE = open_summary_file(*out_folder, "nanomux.log");
    if (LOG_FILE == NULL) goto cleanup;

    S_FILE = open_summary_file(*out_folder, "nanomux_matches.csv");
    if (S_FILE == NULL) goto cleanup;

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

    if (ferror(LOG_FILE)) {
        nob_log(NOB_ERROR, "Failed to write nanomux log: %s", strerror(errno));
        goto cleanup;
    }
    
    if (fprintf(S_FILE, "barcode,matches\n") < 0) {
        nob_log(NOB_ERROR, "Failed to write barcode matches summary: %s", strerror(errno));
        goto cleanup;
    }

    for (size_t i = 0; i < barcodes.count; i++) {
        size_t bc_count = barcodes.items[i].counter;
        char *bc_name = barcodes.items[i].name;
        if (fprintf(S_FILE, "%s,%zu\n", bc_name, bc_count) < 0) {
            nob_log(NOB_ERROR, "Failed to write barcode matches summary: %s", strerror(errno));
            goto cleanup;
        }
        printf("%s: %zu\n", bc_name, bc_count);

        // remove file if empty
        if (bc_count == 0) {
            char *bc_file = barcodes.items[i].out_name;
            if (nob_file_exists(bc_file) && !nob_delete_file(bc_file)) goto cleanup;
        }
    }

    if (fflush(LOG_FILE) != 0 || fflush(S_FILE) != 0) {
        nob_log(NOB_ERROR, "Failed to flush summary files: %s", strerror(errno));
        goto cleanup;
    }

    result = 0;

cleanup:
    // ----------------- CLEAN-UP ---------------------------
    if (thpool != NULL) thpool_destroy(thpool);

    for (size_t i = 0; i < reads.count; i++) free_read(reads.items[i]);
    nob_da_free(reads);

    if (seq != NULL) kseq_destroy(seq);
    if (fp != NULL) {
        int close_status = gzclose(fp);
        if (close_status != Z_OK) {
            nob_log(NOB_ERROR, "Failed to close FASTQ input: %s", zError(close_status));
            result = 1;
        }
    }

    if (S_FILE != NULL && fclose(S_FILE) != 0) {
        nob_log(NOB_ERROR, "Failed to close barcode matches summary: %s", strerror(errno));
        result = 1;
    }
    if (LOG_FILE != NULL && fclose(LOG_FILE) != 0) {
        nob_log(NOB_ERROR, "Failed to close nanomux log: %s", strerror(errno));
        result = 1;
    }

    for (size_t i = 0; i < barcodes.count; i++) free_barcode(&barcodes.items[i]);
    nob_da_free(barcodes);
    nob_sb_free(sb);

    if (output_gate_initialized && !output_gate_destroy(&output_gate)) {
        result = 1;
    }

    if (error_mutex_initialized) {
        int destroy_status = pthread_mutex_destroy(&error_mutex);
        if (destroy_status != 0) {
            nob_log(NOB_ERROR, "Failed to destroy worker error mutex: %s", strerror(destroy_status));
            result = 1;
        }
    }

    if (result == 0) {
        printf("\n");
        nob_log(NOB_INFO, "nanomux done!\n");
    }

    return result;
}
