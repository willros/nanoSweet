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
    size_t read_index;
    size_t score;
    int start;
    int end;
} Candidate_Match;

typedef struct {
    Candidate_Match *items;
    size_t count;
    size_t capacity;
} Candidate_Matches;

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
    Candidate_Matches *candidates;
    size_t barcode_pos;
    size_t k;
    bool trim;
    int barcode_schema;
    pthread_mutex_t *error_mutex;
    bool *worker_failed;
} Match_Thread_Data;

typedef struct {
    Barcode *barcode;
    Pending_Writes *writes;
    pthread_mutex_t *error_mutex;
    bool *worker_failed;
    Output_Gate *output_gate;
} Write_Thread_Data;

typedef struct {
    bool found;
    bool tied;
    size_t barcode_index;
    size_t score;
    int start;
    int end;
} Read_Assignment;

typedef struct {
    size_t assigned;
    size_t ambiguous;
    size_t unclassified;
} Assignment_Stats;

static void mark_worker_failed(pthread_mutex_t *error_mutex, bool *worker_failed)
{
    pthread_mutex_lock(error_mutex);
    *worker_failed = true;
    pthread_mutex_unlock(error_mutex);
}

static bool has_worker_failed(pthread_mutex_t *error_mutex, bool *worker_failed)
{
    pthread_mutex_lock(error_mutex);
    bool failed = *worker_failed;
    pthread_mutex_unlock(error_mutex);
    return failed;
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

static void queue_candidate(Candidate_Matches *matches, size_t read_index,
                            size_t score, int start, int end)
{
    Candidate_Match match = {
        .read_index = read_index,
        .score = score,
        .start = start,
        .end = end,
    };
    nob_da_append(matches, match);
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

static bool write_pending_reads(Write_Thread_Data *td)
{
    if (td->writes->count == 0) return true;

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

    for (size_t i = 0; i < td->writes->count; i++) {
        Pending_Write *pending = &td->writes->items[i];
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

static void process_pending_writes(void *arg)
{
    Write_Thread_Data *td = (Write_Thread_Data *)arg;
    if (!write_pending_reads(td)) {
        mark_worker_failed(td->error_mutex, td->worker_failed);
    }
    free(td);
}

static void process_barcode_matches(void *arg)
{
    Match_Thread_Data *td = (Match_Thread_Data *)arg;
    Barcode *barcode = td->barcode;

    if (!(td->barcode_schema == 1 || td->barcode_schema == 2)) {
        nob_log(NOB_ERROR, "Wrong barcode schema");
        mark_worker_failed(td->error_mutex, td->worker_failed);
        free(td);
        return;
    }

    for (size_t i = 0; i < td->reads->count; i++) {
        Read *read = &td->reads->items[i];
        bool found = false;
        size_t score = 0;
        int start = 0;
        int end = (int)read->len;

        if (td->barcode_schema == 1) {
            Levenshtein_Match forward_match = {0};
            Levenshtein_Match reverse_match = {0};
            bool has_forward = find_best_levenshtein_match(
                read->first_slice, td->barcode_pos,
                barcode->fw, barcode->fw_length, td->k, &forward_match);
            bool has_reverse = false;
            if (!has_forward || forward_match.distance > 0) {
                has_reverse = find_best_levenshtein_match(
                    read->last_slice, td->barcode_pos,
                    barcode->fw_comp, barcode->fw_length, td->k, &reverse_match);
            }

            if (has_forward && (!has_reverse || forward_match.distance <= reverse_match.distance)) {
                found = true;
                score = forward_match.distance;
                if (td->trim) start = forward_match.end;
            } else if (has_reverse) {
                int slice_end = (int)read->len - (int)td->barcode_pos
                              + reverse_match.end - (int)barcode->fw_length;
                if (slice_end > 0) {
                    found = true;
                    score = reverse_match.distance;
                    if (td->trim) end = slice_end;
                }
            }
        } else {
            Levenshtein_Match left_forward = {0};
            Levenshtein_Match right_reverse = {0};
            Levenshtein_Match left_reverse = {0};
            Levenshtein_Match right_forward = {0};

            bool has_forward_pair = find_best_levenshtein_match(
                read->first_slice, td->barcode_pos,
                barcode->fw, barcode->fw_length, td->k, &left_forward);
            if (has_forward_pair) {
                has_forward_pair = find_best_levenshtein_match(
                    read->last_slice, td->barcode_pos,
                    barcode->rv_comp, barcode->rv_length, td->k, &right_reverse);
            }
            size_t forward_score = has_forward_pair
                ? left_forward.distance + right_reverse.distance : 0;

            bool has_reverse_pair = false;
            if (!has_forward_pair || forward_score > 0) {
                has_reverse_pair = find_best_levenshtein_match(
                    read->first_slice, td->barcode_pos,
                    barcode->rv, barcode->rv_length, td->k, &left_reverse);
                if (has_reverse_pair) {
                    has_reverse_pair = find_best_levenshtein_match(
                        read->last_slice, td->barcode_pos,
                        barcode->fw_comp, barcode->fw_length, td->k, &right_forward);
                }
            }

            size_t reverse_score = has_reverse_pair
                ? left_reverse.distance + right_forward.distance : 0;

            if (has_forward_pair && (!has_reverse_pair || forward_score <= reverse_score)) {
                int slice_end = (int)read->len - (int)td->barcode_pos
                              + right_reverse.end - (int)barcode->rv_length;
                if (slice_end > 0) {
                    found = true;
                    score = forward_score;
                    if (td->trim) {
                        start = left_forward.end;
                        end = slice_end;
                    }
                }
            } else if (has_reverse_pair) {
                int slice_end = (int)read->len - (int)td->barcode_pos
                              + right_forward.end - (int)barcode->fw_length;
                if (slice_end > 0) {
                    found = true;
                    score = reverse_score;
                    if (td->trim) {
                        start = left_reverse.end;
                        end = slice_end;
                    }
                }
            }
        }

        if (found && start < end) {
            queue_candidate(td->candidates, i, score, start, end);
        }
    }

    free(td);
}

static bool resolve_candidates(Barcodes *barcodes, Reads *reads,
                               Candidate_Matches *candidate_sets,
                               Pending_Writes *writes,
                               Assignment_Stats *stats)
{
    Read_Assignment *assignments = calloc(reads->count, sizeof(*assignments));
    if (assignments == NULL && reads->count > 0) {
        nob_log(NOB_ERROR, "Failed to allocate read assignments");
        return false;
    }

    for (size_t barcode_index = 0; barcode_index < barcodes->count; barcode_index++) {
        Candidate_Matches *matches = &candidate_sets[barcode_index];
        for (size_t i = 0; i < matches->count; i++) {
            Candidate_Match *candidate = &matches->items[i];
            Read_Assignment *assignment = &assignments[candidate->read_index];

            if (!assignment->found || candidate->score < assignment->score) {
                assignment->found = true;
                assignment->tied = false;
                assignment->barcode_index = barcode_index;
                assignment->score = candidate->score;
                assignment->start = candidate->start;
                assignment->end = candidate->end;
            } else if (candidate->score == assignment->score &&
                       barcode_index != assignment->barcode_index) {
                assignment->tied = true;
            }
        }
    }

    for (size_t read_index = 0; read_index < reads->count; read_index++) {
        Read_Assignment *assignment = &assignments[read_index];
        if (!assignment->found) {
            stats->unclassified++;
        } else if (assignment->tied) {
            stats->ambiguous++;
        } else {
            Barcode *barcode = &barcodes->items[assignment->barcode_index];
            queue_barcode_read(&writes[assignment->barcode_index],
                               &reads->items[read_index],
                               assignment->start, assignment->end);
            barcode->counter++;
            stats->assigned++;
        }
    }

    free(assignments);
    return true;
}

static bool dispatch_barcodes(threadpool thpool, Barcodes *barcodes, Reads *reads,
                              size_t barcode_pos, size_t k, bool trim,
                              int barcode_schema, pthread_mutex_t *error_mutex,
                              bool *worker_failed, Output_Gate *output_gate,
                              Assignment_Stats *stats)
{
    if (barcodes->count == 0) {
        stats->unclassified += reads->count;
        return true;
    }

    bool success = false;
    bool all_jobs_scheduled = true;
    Candidate_Matches *candidate_sets = calloc(barcodes->count, sizeof(*candidate_sets));
    Pending_Writes *writes = calloc(barcodes->count, sizeof(*writes));
    if (candidate_sets == NULL || writes == NULL) {
        nob_log(NOB_ERROR, "Failed to allocate batch match data");
        goto cleanup;
    }

    for (size_t i = 0; i < barcodes->count; i++) {
        Match_Thread_Data *td = malloc(sizeof(*td));
        if (td == NULL) {
            nob_log(NOB_ERROR, "Failed to allocate match thread data");
            all_jobs_scheduled = false;
            break;
        }

        td->barcode = &barcodes->items[i];
        td->reads = reads;
        td->candidates = &candidate_sets[i];
        td->barcode_pos = barcode_pos;
        td->k = k;
        td->trim = trim;
        td->barcode_schema = barcode_schema;
        td->error_mutex = error_mutex;
        td->worker_failed = worker_failed;

        if (thpool_add_work(thpool, process_barcode_matches, td) != 0) {
            nob_log(NOB_ERROR, "Failed to schedule barcode matching");
            free(td);
            all_jobs_scheduled = false;
            break;
        }
    }

    thpool_wait(thpool);
    if (!all_jobs_scheduled || has_worker_failed(error_mutex, worker_failed)) {
        goto cleanup;
    }

    if (!resolve_candidates(barcodes, reads, candidate_sets, writes, stats)) {
        goto cleanup;
    }

    for (size_t i = 0; i < barcodes->count; i++) {
        if (writes[i].count == 0) continue;

        Write_Thread_Data *td = malloc(sizeof(*td));
        if (td == NULL) {
            nob_log(NOB_ERROR, "Failed to allocate write thread data");
            all_jobs_scheduled = false;
            break;
        }

        td->barcode = &barcodes->items[i];
        td->writes = &writes[i];
        td->error_mutex = error_mutex;
        td->worker_failed = worker_failed;
        td->output_gate = output_gate;

        if (thpool_add_work(thpool, process_pending_writes, td) != 0) {
            nob_log(NOB_ERROR, "Failed to schedule barcode output");
            free(td);
            all_jobs_scheduled = false;
            break;
        }
    }

    thpool_wait(thpool);
    if (!all_jobs_scheduled || has_worker_failed(error_mutex, worker_failed)) {
        goto cleanup;
    }

    success = true;

cleanup:
    if (candidate_sets != NULL) {
        for (size_t i = 0; i < barcodes->count; i++) {
            nob_da_free(candidate_sets[i]);
        }
    }
    if (writes != NULL) {
        for (size_t i = 0; i < barcodes->count; i++) {
            nob_da_free(writes[i]);
        }
    }
    free(candidate_sets);
    free(writes);
    return success;
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
    Assignment_Stats assignment_stats = {0};
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
                                   *trim, barcode_schema, &error_mutex, &worker_failed, &output_gate, &assignment_stats)) {
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
                               *trim, barcode_schema, &error_mutex, &worker_failed, &output_gate, &assignment_stats)) {
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
    printf("INFO: Assigned reads: %zu reads\n", assignment_stats.assigned);
    printf("INFO: Ambiguous reads: %zu reads\n", assignment_stats.ambiguous);
    printf("INFO: Unclassified reads: %zu reads\n", assignment_stats.unclassified);
    fprintf(LOG_FILE, "Processed %zu reads\n", counter);
    fprintf(LOG_FILE, "Reads shorter than p: %zu reads\n", reads_shorter_than_p);
    fprintf(LOG_FILE, "Assigned reads: %zu reads\n", assignment_stats.assigned);
    fprintf(LOG_FILE, "Ambiguous reads: %zu reads\n", assignment_stats.ambiguous);
    fprintf(LOG_FILE, "Unclassified reads: %zu reads\n", assignment_stats.unclassified);

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
