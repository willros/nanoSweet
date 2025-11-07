#define COMMON_IMPLEMENTATION
#include "common.h"
#define FLAG_IMPLEMENTATION
#include "flag.h"
#include <zlib.h>
#include "kseq.h"
#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include "thpool.h"
#include <string.h>

KSEQ_INIT(gzFile, gzread)
#define READ_BUFFER (2 * 1000)

typedef struct {
    size_t min_qual;
    size_t min_len;
    size_t max_len;

    const char *in_file;
    const char *out_file;

    size_t raw_reads;
    size_t too_short;
    size_t too_long;
    size_t too_bad;
    size_t qualified_reads;
} Fastq_File; 

typedef struct {
    Fastq_File *items;
    size_t count;
    size_t capacity;
} Fastq_Files;

typedef struct {
    Fastq_File *f;
    Reads *reads;
    size_t start;
    size_t end;
    gzFile out_file;
    pthread_mutex_t *print_mutex;
} Thread_Data;


bool parse_input( const char *input, const char *output, Fastq_Files *fastq_files, size_t min_qual, size_t min_len, size_t max_len ) {
    Nob_File_Type type = nob_get_file_type(input);
    Nob_File_Paths files = {0};

    switch (type) {
        case NOB_FILE_DIRECTORY: {
            nob_log(NOB_INFO, "%s is a directory", input);

            if (!nob_read_entire_dir(input, &files)) {
                nob_log(NOB_ERROR, "Failed to read directory %s", input);
                return false;
            }

            char real_path[512];
            char outfile[512];

            for (size_t i = 0; i < files.count; ++i) {
                const char *file = files.items[i];
                if (strcmp(file, ".") == 0) continue;
                if (strcmp(file, "..") == 0) continue;
                if (*file == '.') continue;
                if (!is_fastq(file)) continue;

                // realpath and outfile
                snprintf(real_path, sizeof(real_path), "%s/%s", input, file);
                snprintf(outfile, sizeof(outfile), "%s/%s_nanotrim.fq.gz", output, file);

                Fastq_File fastq_file = {
                    .min_qual = min_qual,
                    .min_len = min_len,
                    .max_len = max_len,
                    .in_file = strdup(real_path),
                    .out_file = strdup(outfile),
                    .raw_reads = 0,
                    .too_short = 0,
                    .too_long = 0,
                    .too_bad = 0,
                    .qualified_reads = 0,
                };
                nob_da_append(fastq_files, fastq_file);
            }

            break;
        }

        case NOB_FILE_REGULAR: {
            if (!is_fastq(input)) {
                nob_log(NOB_ERROR, "%s is not a fastq file.", input);
                return false;
            }

            nob_log(NOB_INFO, "%s is a file", input);

            char outfile[512];
            char *base_name = basename(input);
            snprintf(outfile, sizeof(outfile), "%s/%s.filtered", output, base_name);

            Fastq_File fastq_file = {
                .min_qual = min_qual,
                .min_len = min_len,
                .max_len = max_len,
                .in_file = strdup(input),
                .out_file = strdup(outfile),
                .raw_reads = 0,
                .too_short = 0,
                .too_long = 0,
                .too_bad = 0,
                .qualified_reads = 0,
            };
            nob_da_append(fastq_files, fastq_file);

            break;
        }

        default:
            nob_log(NOB_ERROR, "input: %s has an unknown type", input);
            return false;
    }

    nob_da_free(files);
    return true;
}


void parse_fastq(void *arg) 
{
    Thread_Data *td = (Thread_Data *)arg;
    size_t local_raw = 0;
    size_t local_short = 0;
    size_t local_long = 0;
    size_t local_bad = 0;
    size_t local_passed = 0;

    Reads *reads = td->reads;
    Fastq_File *f = td->f;

    for (size_t idx = td->start; idx < td->end && idx < reads->count; idx++) {
        Read cur_read = reads->items[idx];
        local_raw++;

        if (cur_read.len < f->min_len) { local_short++; continue; }
        if (cur_read.len > f->max_len) { local_long++; continue; }

        if (average_qual(cur_read.qual, cur_read.len) < (double)f->min_qual) { local_bad++; continue; }

        pthread_mutex_lock(td->print_mutex);
            if (!append_read_to_gzip_fastq(td->out_file, &cur_read, 0, cur_read.len)) {
                // If append fails, unlock then exit
                pthread_mutex_unlock(td->print_mutex);
                exit(1);
            }
            f->qualified_reads++;
        pthread_mutex_unlock(td->print_mutex);
        local_passed++;
    }

    pthread_mutex_lock(td->print_mutex);
        f->raw_reads += local_raw;
        f->too_short += local_short;
        f->too_long += local_long;
        f->too_bad += local_bad;
    pthread_mutex_unlock(td->print_mutex);

    free(td);
}


int main(int argc, char **argv) {

    // flag.h arguments
    char **input = flag_str("f", "", "Path to input folder or file (MANDATORY)");
    char **out_dir = flag_str("o", "", "Name of output folder (MANDATORY)");
    size_t *min_len = flag_size("r", 0, "Minimum read length");
    size_t *max_len = flag_size("R", 1000*1000, "Maximum read length");
    size_t *min_qual = flag_size("q", 0, "Minimum quality");
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

    if (strcmp(*input, "") == 0 || strcmp(*out_dir, "") == 0) {
        nob_log(NOB_ERROR, "At least one of the mandatory arguments are missing");
        flag_print_options(stderr);
        return 1;
    }
    
    if (!nob_mkdir_if_not_exists(*out_dir)) {
        nob_log(NOB_ERROR, "exiting");
        return 1;
    }

    nob_log(NOB_INFO, "Input:               %20s", *input);
    nob_log(NOB_INFO, "Output:              %20s", *out_dir);
    nob_log(NOB_INFO, "Minimum read length: %20zu", *min_len);
    nob_log(NOB_INFO, "Maximum read length: %20zu", *max_len);
    nob_log(NOB_INFO, "Minimum quality:     %20zu", *min_qual);
    nob_log(NOB_INFO, "Number of threads:   %20zu", *num_threads);


    // -------------- PARSE INPUT ---------------------
    Fastq_Files fastq_files = {0};
    if(!parse_input(*input, *out_dir, &fastq_files, *min_qual, *min_len, *max_len)) return 1;

    // -------------- GENERATE THREAD POOL ---------------------
    nob_log(NOB_INFO, "Generating threadpool with %zu threads", *num_threads);
    threadpool thpool = thpool_init(*num_threads);

    
    // -------------- LOOP THROUGH EVERY INPUT FILE ---------------------
    pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
    Reads reads = {0};

    for (size_t fi = 0; fi < fastq_files.count; fi++) {
        Fastq_File *f = &fastq_files.items[fi];
        gzFile in_file = gzopen(f->in_file, "r"); 
        if (!in_file) {
            nob_log(NOB_ERROR, "Failed to open %s file, exiting", f->in_file);
            return 1;
        }
        gzFile out_file = gzopen(f->out_file, "ab"); 
        if (!out_file) {
            nob_log(NOB_ERROR, "Failed to open %s file, exiting", f->out_file);
            gzclose(in_file);
            return 1;
        }
        kseq_t *seq = kseq_init(in_file); 
        if (seq == NULL) {
            nob_log(NOB_ERROR, "Could not initialize %s file, exiting", f->in_file);
            gzclose(in_file);
            gzclose(out_file);
            return 1;
        }

        // read loop
        while (kseq_read(seq) >= 0) { 
            Read read = {0};
            read.seq = strdup(seq->seq.s);
            read.name = strdup(seq->name.s);
            read.qual = strdup(seq->qual.s);
            read.len = (size_t)seq->seq.l;
            read.first_slice = NULL;
            read.last_slice = NULL;
            
            nob_da_append(&reads, read);
        
            if (reads.count >= READ_BUFFER) {
                size_t total = reads.count;
                size_t per_thread = total / *num_threads;
                size_t rest = total % *num_threads;
                size_t start = 0;
                size_t end = 0;

                for (size_t t = 0; t < *num_threads; t++) {
                    Thread_Data *td = malloc(sizeof(Thread_Data));
                    if (!td) {
                        nob_log(NOB_ERROR, "Failed to allocate thread data");
                        return 1;
                    }
                    
                    start = end;
                    end += per_thread;
                    if (t == (*num_threads - 1)) end += rest; 

                    td->f = f;
                    td->reads = &reads;
                    td->start = start;
                    td->end = end;
                    td->out_file = out_file;
                    td->print_mutex = &print_mutex;
                    
            		thpool_add_work(thpool, parse_fastq, (void *)td);
                }
            	thpool_wait(thpool);

                // Clean up reads
                for (size_t ri = 0; ri < reads.count; ri++) free_read(reads.items[ri]);
                reads.count = 0;
            }
        }

        // ------------- IF ANY READS LEFT -------------------
        if (reads.count > 0) {
            size_t total = reads.count;
            size_t per_thread = total / *num_threads;
            size_t rest = total % *num_threads;
            size_t start = 0;
            size_t end = 0;
            for (size_t t = 0; t < *num_threads; t++) {
                Thread_Data *td = malloc(sizeof(Thread_Data));
                if (!td) {
                    nob_log(NOB_ERROR, "Failed to allocate thread data");
                    return 1;
                }
                
                start = end;
                end += per_thread;
                if (t == (*num_threads - 1)) end += rest;
                td->f = f;
                td->reads = &reads;
                td->start = start;
                td->end = end;
                td->out_file = out_file;
                td->print_mutex = &print_mutex;
                
                thpool_add_work(thpool, parse_fastq, (void *)td);
            }
            thpool_wait(thpool);

            // Clean up reads
            for (size_t ri = 0; ri < reads.count; ri++) free_read(reads.items[ri]);
            reads.count = 0;
        }
        
        kseq_destroy(seq); 
        gzclose(in_file); 
        gzclose(out_file); 
    }


    // -------------- PRINT TO SUMMARY FILES ---------------------
    FILE *LOG_FILE = open_summary_file(*out_dir, "nanotrim_log.csv");
    fprintf(LOG_FILE, "file,raw_reads,passed_reads,short,long,bad_quality\n");
    for (size_t i = 0; i < fastq_files.count; i++) {
        Fastq_File f = fastq_files.items[i];
        nob_log(
            NOB_INFO, 
            "%-10s: %zu raw reads (%zu passed) --> Too short: %-5zu | Too long: %-5zu | Too low quality: %-5zu", 
            f.in_file, f.raw_reads, f.qualified_reads, f.too_short, f.too_long, f.too_bad);
        fprintf(LOG_FILE, "%s,%zu,%zu,%zu,%zu,%zu\n", f.in_file, f.raw_reads, f.qualified_reads, f.too_short, f.too_long, f.too_bad);
    }
    
    // -------------- CLEAN UP ---------------------
	thpool_destroy(thpool);
    pthread_mutex_destroy(&print_mutex);
    nob_da_free(fastq_files);
    fclose(LOG_FILE);

    return 0;
}
