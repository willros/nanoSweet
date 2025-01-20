#define NOB_IMPLEMENTATION
#include "nob.h"
#include "kseq.h"
#include <stdio.h>
#include <zlib.h>
#include <pthread.h>
#include "thpool.h"
#include <stdint.h>


uint32_t djb2(const char *buf, size_t buf_size) {
    uint32_t hash = 5381;

    for (size_t i = 0; i < buf_size; ++i) {
        hash = ((hash << 5) + hash) + (uint32_t)buf[i]; 
    }

    return hash;
}


char *basename(char const *path) {
    char *s = strrchr(path, '/');
    if (!s)
        return strdup(path);
    else
        return strdup(s + 1);
}


#define hash_init(ht, cap) \
    do { \
        (ht)->items = malloc(sizeof(*(ht)->items)*cap); \
        memset((ht)->items, 0, sizeof(*(ht)->items)*cap); \
        (ht)->count = 0; \
        (ht)->capacity = (cap); \
    } while(0)

#define hash_resize(ht) \
    do { \
        size_t new_capacity = (ht)->capacity * 2; \
        Read *new_items = malloc(sizeof(*(ht)->items) * new_capacity); \
        if (!new_items) { \
            nob_log(NOB_ERROR, "Failed to allocate memory for hash table resize"); \
            break; \
        } \
        memset(new_items, 0, sizeof(*(ht)->items) * new_capacity); \
        for (size_t i = 0; i < (ht)->capacity; i++) { \
            if ((ht)->items[i].occupied) { \
                uint32_t h = djb2((ht)->items[i].key, strlen((ht)->items[i].key)) % new_capacity; \
                while (new_items[h].occupied) { \
                    h = (h + 1) % new_capacity; \
                } \
                new_items[h] = (ht)->items[i]; \
            } \
        } \
        free((ht)->items); \
        (ht)->items = new_items; \
        (ht)->capacity = new_capacity; \
    } while(0)

typedef struct {
    const char *key;
    int value;
    bool occupied;
} Read;

typedef struct {
    Read *items;
    size_t count;
    size_t capacity;
} Reads;


bool is_fastq(const char *file) {
    return strstr(file, "fastq") || strstr(file, "fq") || strstr(file, "fa");
}

bool must_be_digit(const char *arg) {
    for(; *arg != '\0'; arg++) {
        if (!isdigit(*arg)) return false;
    }
    return true;
}

int append_read_to_gzip_fastq(gzFile gzfp, const char *name, const char *seq, const char *qual) {
    int result = 0;
    size_t max_len = strlen(seq);
    size_t buffer_size = max_len + 10; 
    char *buffer = (char *)malloc(buffer_size);

    if (!buffer) {
        nob_log(NOB_ERROR, "Memory allocation failed");
        return 1;
    }

    int len = snprintf(buffer, buffer_size, "@%s\n", name);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%s\n", seq);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write to gzip file");
        nob_return_defer(1);
    }

    if (qual) {
        len = snprintf(buffer, buffer_size, "+\n");
        if (gzwrite(gzfp, buffer, len) != len) {
            nob_log(NOB_ERROR, "Failed to write to gzip file");
            nob_return_defer(1);
        }

        len = snprintf(buffer, buffer_size, "%s\n", qual);
            if (gzwrite(gzfp, buffer, len) != len) {
                nob_log(NOB_ERROR, "Failed to write to gzip file");
                nob_return_defer(1);
            }
    }

defer:
    free(buffer);
    return result;
}

typedef struct {
    char *in_file;
    char *clean_file;
    char *log_file_file;
    FILE *log_file_all;
    pthread_mutex_t *log_file_all_mutex;
} File; 

typedef struct {
    File *items;
    size_t count;
    size_t capacity;
} Files;



KSEQ_INIT(gzFile, gzread)
bool nanodup_file(File *file) {
    gzFile fp = gzopen(file->in_file, "r"); 
    if (!fp) {
        nob_log(NOB_INFO, "Failed to open fastq file: '%s'", file->in_file);
        return false;
    }
    gzFile out_file = gzopen(file->clean_file, "wb"); 
    if (!out_file) {
        nob_log(NOB_ERROR, "Failed to open %s file, exiting", file->clean_file);
        return false;
    }

    Reads ht = {0};
    hash_init(&ht, 1024*10);

    int l;
    kseq_t *seq = kseq_init(fp); 
    int num_reads = 0;

    while ((l = kseq_read(seq)) >= 0) { 

        if (ht.count >= (ht.capacity * 0.7)) { 
            hash_resize(&ht);
        }

        uint32_t h = djb2(seq->seq.s, strlen(seq->seq.s))%ht.capacity;
        for (size_t i = 0; i < ht.capacity && ht.items[h].occupied && strcmp(ht.items[h].key, seq->seq.s) != 0; ++i) {
            h = (h + 1)%ht.capacity;
        }

        if (ht.items[h].occupied) {
            if (strcmp(ht.items[h].key, seq->seq.s) != 0) {
                nob_log(NOB_ERROR, "Table overflow");
                return false;
            }
            ht.items[h].value += 1;
        } else {
            ht.items[h].occupied = true;
            ht.items[h].key = strdup(seq->seq.s);
            ht.items[h].value = 1;
            ht.count++;
            append_read_to_gzip_fastq(out_file, seq->name.s, seq->seq.s, seq->qual.s);
        }
        num_reads += 1;
    }

    Reads freq = {0};
    for (size_t i = 0; i < ht.capacity; ++i) {
        if (ht.items[i].occupied) {
            nob_da_append(&freq, ht.items[i]);
        }
    }

    FILE *log_file_file = fopen(file->log_file_file, "ab");
    fprintf(log_file_file, "read,count\n");

    for (size_t i = 0; i < freq.count; ++i) {
        if (freq.items[i].value > 1) {
            fprintf(log_file_file, "%s,%i\n", freq.items[i].key, freq.items[i].value);
        }
    }

    // log file all
    pthread_mutex_lock(file->log_file_all_mutex);
    fprintf(file->log_file_all, "%s,%i,%i,%i\n", file->in_file, num_reads, (int)freq.count, num_reads - (int)freq.count);
    pthread_mutex_unlock(file->log_file_all_mutex);

    nob_log(NOB_INFO, "%s contained: %i duplicates", file->in_file, num_reads - (int)freq.count);

    kseq_destroy(seq); 
    gzclose(fp); 
    gzclose(out_file);
    fclose(log_file_file);
    nob_da_free(ht);

    return true;
}

void task_parse_fastq_file(void *arg) {
    File *file = (File *)arg;

    if (!nanodup_file(file)) {
        return;
    }
}

char *usage = 
"[USAGE]: nanodup -i <input> -o <output> [options]\n"
"   -i    <input>             Path of folder or file\n"
"   -o    <output>            Name of output folder.\n"
"   -t    <threads>           Number of threads to use. Optional: Default 1\n";

int main(int argc, char **argv) {

    // required args
    bool i_arg = false;
    char *input;

    char *output;
    bool o_arg = false;

    char c;
    int t_arg = 1;
    while ((c = getopt(argc, argv, "i:o:t:")) != -1) {
        switch (c) {
            case 'i':
                i_arg = true;
                input = optarg;
                break;
            case 'o':
                o_arg = true;
                output = optarg;
                break;
            case 't':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-t must be digit");
                    nob_log(NOB_ERROR, "%s", usage);
                    return 1;
                }
                t_arg = atoi(optarg);
                break;
            }
    }
    if (!(i_arg)) {
        nob_log(NOB_ERROR, "You must provide the input path");
        nob_log(NOB_ERROR, "%s", usage);
        return 1;
    }
    if (!(o_arg)) {
        nob_log(NOB_ERROR, "You must provide the output path");
        nob_log(NOB_ERROR, "%s", usage);
        return 1;
    }

    nob_log(NOB_INFO, "Input:               %20s", input);
    nob_log(NOB_INFO, "Output:              %20s", output);
    nob_log(NOB_INFO, "Number of threads:   %20i", t_arg);

    if (!nob_mkdir_if_not_exists(output)) {
        nob_log(NOB_ERROR, "exiting");
        return 1;
    }

    Nob_File_Type type = nob_get_file_type(input);
    Nob_File_Paths files = {0};
    Files fastq_files = {0};

    char log_file_all[256]; 
    snprintf(log_file_all, sizeof(log_file_all), "%s/nanodup_log.csv", output);
    FILE *LOG_FILE_ALL = fopen(log_file_all, "ab");
    if (LOG_FILE_ALL == NULL) {
        nob_log(NOB_ERROR, "Could create log file");
        return 1;
    }
    fprintf(LOG_FILE_ALL, "file,num_raw,num_unique,num_duplicated\n");
    pthread_mutex_t log_file_mutex = PTHREAD_MUTEX_INITIALIZER;

    switch (type) {
        case NOB_FILE_DIRECTORY: {
            nob_log(NOB_INFO, "`%s` is a directory", input);
            if (!nob_read_entire_dir(input, &files)) {
                nob_log(NOB_ERROR, "Failed to read directory `%s`", input);
                return 1;
            }


            char clean_file[512];
            char file_log_file[512];
            for (size_t i = 0; i < files.count; ++i) {
                const char *file = files.items[i];

                if (strcmp(file, ".") == 0) continue;
                if (strcmp(file, "..") == 0) continue;
                if (*file == '.') continue;
                if (!is_fastq(file)) continue;

                char *base_name = basename(file);
                char in_file[1024];
                snprintf(in_file, sizeof(in_file), "%s/%s", input, file);

                snprintf(clean_file, sizeof(clean_file), "%s/%s.nanoduped.fq.gz", output, base_name);
                snprintf(file_log_file, sizeof(file_log_file), "%s/%s.nanodup.log", output, base_name);

                File fastq_file = { 
                    .in_file = strdup(in_file),
                    .clean_file = strdup(clean_file),
                    .log_file_file = strdup(file_log_file),
                    .log_file_all = LOG_FILE_ALL,
                    .log_file_all_mutex = &log_file_mutex,
                };
                nob_da_append(&fastq_files, fastq_file);
            }
            break;
        }
        case NOB_FILE_REGULAR: {
            if (!is_fastq(input)) {
                nob_log(NOB_ERROR, "`%s` is not a fastq file...", input);
                return 1;
            }
            nob_log(NOB_INFO, "`%s` is a file", input);
            char *base_name = basename(input);

            char clean_file[1024];
            snprintf(clean_file, sizeof(clean_file), "%s/%s.nanoduped.fq.gz", output, base_name);

            char file_log_file[1024];
            snprintf(file_log_file, sizeof(file_log_file), "%s/%s.duplicated", output, base_name);

            File fastq_file = { 
                .in_file = strdup(input),
                .clean_file = strdup(clean_file),
                .log_file_file = strdup(file_log_file),
                .log_file_all = LOG_FILE_ALL,
                .log_file_all_mutex = &log_file_mutex,
            };

            nob_da_append(&fastq_files, fastq_file);
            break;
        }
        default:
            nob_log(NOB_ERROR, "input: `%s` has an unknown type", input);
			return 1;
    }

    nob_log(NOB_INFO, "Generating threadpool with %i threads", t_arg);
    threadpool thpool = thpool_init(t_arg);

    for (int i = 0; i < fastq_files.count; i++) {
		thpool_add_work(thpool, task_parse_fastq_file, (void *)&fastq_files.items[i]);
    }

	thpool_wait(thpool);
	thpool_destroy(thpool);
    pthread_mutex_destroy(&log_file_mutex);
    nob_log(NOB_INFO, "nanodup done!");

    nob_da_free(fastq_files);
    nob_da_free(files);

    fclose(LOG_FILE_ALL);

    return 0;

}
