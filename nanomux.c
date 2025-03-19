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
#include <zlib.h>
#include "kseq.h"
#include <limits.h> 
#include <stdint.h>
#include <pthread.h>
#include "thpool.h"


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

    int end_pos = -1;
    for (int j = needle_len; j <= haystack_len; j++) {
        if (dp[needle_len][j] <= k) {
            return j;
        }
    }
    return end_pos;
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
    int length;
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


// TODO: update the adapter len when changing the adapter
// TODO: change the parameters of the adapter position and so on
const char *ADAPTER = "TACTTCGTTCAGTTACGTATTGCT";
#define ADAPTER_LEN 24
#define CONCATENATED_READ_LEN_MIN 2500
#define CONCATENATED_READ_LEN_MAX 5000
#define ADAPTER_POSITION 1200

void split_reads_by_adapter(const Read *read, Reads *reads, int *counter, int read_len_min, int read_len_max) {

    int match = levenshtein_distance(read->seq, ADAPTER, 1); // k = 1
    
    if (match > ADAPTER_POSITION) {  
        *counter += 1;
        // two loops, one for each new read
        for (int i = 1; i < 3; i++) {
            
            int start = (i == 1) ? 0 : match;
            int end = (i == 1) ? match : read->length;
            int new_length = end - start;
            
            if (new_length >= read_len_min && new_length <= read_len_max) {
                char new_read_name[100];
                sprintf(new_read_name, "%s_%i", read->name, i);

                char new_read_seq[3500];
                snprintf(new_read_seq, sizeof(new_read_seq), "%.*s", (int)new_length, read->seq + start);

                char new_read_qual[3500];
                snprintf(new_read_qual, sizeof(new_read_qual), "%.*s", (int)new_length, read->qual + start);

                // add everything to the da Reads
                Read new_read = {0};
                new_read.seq = strdup(new_read_seq);
                new_read.name = strdup(new_read_name);
                new_read.qual = strdup(new_read_qual);
                new_read.length = new_length;
                nob_da_append(reads, new_read);
            }
        }
    }
}

KSEQ_INIT(gzFile, gzread)
Reads parse_fastq(const char *fastq_file_path, int read_len_min, int read_len_max, bool split_reads) {
    
    gzFile fp = gzopen(fastq_file_path, "r"); 
    if (!fp) {
        nob_log(NOB_INFO, "Failed to open fastq file, exiting");
        exit(1);
    }
    
    kseq_t *seq = kseq_init(fp); 
    int l;
    Reads reads = {0};
    
    int counter = 0;
    
    while ((l = kseq_read(seq)) >= 0) { 
        int length = strlen(seq->seq.s);
        
        if (split_reads) {
            if (length >= CONCATENATED_READ_LEN_MIN && length <= CONCATENATED_READ_LEN_MAX) {
                Read concated_read = {0};
                concated_read.seq = strdup(seq->seq.s);
                concated_read.name = strdup(seq->name.s);
                concated_read.qual = strdup(seq->qual.s);
                concated_read.length = length;

                split_reads_by_adapter(&concated_read, &reads, &counter, read_len_min, read_len_max);
                free_read(&concated_read);
                continue;
            }
        }
        
        if (length >= read_len_min && length <= read_len_max) {
            Read read = {0};
            read.seq = strdup(seq->seq.s);
            read.name = strdup(seq->name.s);
            read.qual = strdup(seq->qual.s);
            read.length = length;
            nob_da_append(&reads, read);
        }
    }
    
    if (split_reads) nob_log(NOB_INFO, "Number of splitted reads: %i", counter);
    
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

void complement_sequence(const char *seq, char *comp_seq, size_t length) {
    for (size_t i = 0; i < length; i++) {
        comp_seq[length - 1 - i] = complement(seq[i]);
    }
    comp_seq[length] = '\0';
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

void substring(const char *source, size_t start, size_t length, char *dest) {
    memcpy(dest, source + start, length);
    dest[length] = '\0';
}

int append_read_to_gzip_fastq_trim(gzFile gzfp, Read *read, int start, int end) {
    int result = 0;

    if (start < 0) start = 0;
    if (end > (int)read->length) end = read->length;
    size_t trimmed_length = end - start;

    size_t buffer_size = trimmed_length + 10; 
    char *buffer = (char *)malloc(buffer_size);

    if (!buffer) {
        nob_log(NOB_ERROR, "Memory allocation failed");
        return 1;
    }

    int len = snprintf(buffer, buffer_size, "@%s\n", read->name);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write name to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%.*s\n", (int)trimmed_length, read->seq + start);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write trimmed sequence to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "+\n");
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write separator to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%.*s\n", (int)trimmed_length, read->qual + start);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write trimmed quality to gzip file");
        nob_return_defer(1);
    }

defer:
    free(buffer);
    return result;
}

int append_read_to_gzip_fastq(gzFile gzfp, Read *read) {
    int result = 0;
    size_t max_len = strlen(read->seq);
    size_t buffer_size = max_len + 10; 
    char *buffer = (char *)malloc(buffer_size);

    if (!buffer) {
        nob_log(NOB_ERROR, "Memory allocation failed");
        return 1;
    }

    int len = snprintf(buffer, buffer_size, "@%s\n", read->name);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write name to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%s\n", read->seq);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write sequence to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "+\n");
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write separator to gzip file");
        nob_return_defer(1);
    }

    len = snprintf(buffer, buffer_size, "%s\n", read->qual);
    if (gzwrite(gzfp, buffer, len) != len) {
        nob_log(NOB_ERROR, "Failed to write quality to gzip file");
        nob_return_defer(1);
    }

defer:
    free(buffer);
    return result;
}

void process_single_barcode(
    const Barcode b,
    const Reads reads,
    const char *output,
    const size_t barcode_pos,
    const int k,
    const int trim,
    FILE *S_FILE,
    pthread_mutex_t *s_mutex
) {
    int counter = 0;
    
    // Save fastq
    char fastq_name[256]; 
    snprintf(fastq_name, sizeof(fastq_name), "%s/%s.fq.gz", output, b.name);
    
    gzFile new_fastq = gzopen(fastq_name, "ab");
    if (!new_fastq) {
        nob_log(NOB_INFO, "Failed to open new fastq_file for appending, exiting");
        return;
    }
    
    // TODO: change this to read_len_max
    char first_part[1000];
    char last_part[1000];
    
    for (size_t j = 0; j < reads.count; ++j) {
        
        Read r = reads.items[j];
        size_t first_part_start = 0;
        substring(r.seq, first_part_start, barcode_pos, first_part); 
        size_t last_part_start = r.length - barcode_pos;
        substring(r.seq, last_part_start, barcode_pos, last_part); 
        
        int match_first_fw = levenshtein_distance(first_part, b.fw, k);
        // fw ------ revcomp(rv)
        if (match_first_fw != -1) {
            int match_last_fw = levenshtein_distance(last_part, b.rv_comp, k);
            if (match_last_fw != -1) {
                counter++;
                int trim_right = last_part_start + match_last_fw - b.rv_length;
                // trimming or not
                if (trim) {
                    append_read_to_gzip_fastq_trim(new_fastq, &r, match_first_fw, trim_right);
                } else {
                    append_read_to_gzip_fastq(new_fastq, &r);
                }
                continue;
            }
        }
        
        // rv ------ revcomp(fw)
        int match_first_rv = levenshtein_distance(first_part, b.rv, k);
        if (match_first_rv != -1) {
            int match_last_rv = levenshtein_distance(last_part, b.fw_comp, k);
            if (match_last_rv != -1) {
                int trim_right = last_part_start + match_last_rv - b.fw_length;
                counter++;
                // trimming or not
                if (trim) {
                    append_read_to_gzip_fastq_trim(new_fastq, &r, match_first_rv, trim_right);
                } else {
                    append_read_to_gzip_fastq(new_fastq, &r);
                }
            }
        }
    }
    gzclose(new_fastq);
    
    pthread_mutex_lock(s_mutex);
    nob_log(NOB_INFO, "barcode: %30s matches: %i", b.name, counter);
    fprintf(S_FILE, "%s,%i\n", b.name, counter);
    pthread_mutex_unlock(s_mutex);
}


void process_single_barcode_forward(
    const Barcode b,
    const Reads reads,
    const char *output,
    const size_t barcode_pos,
    const int k,
    const int trim,
    FILE *S_FILE,
    pthread_mutex_t *s_mutex
) {
    int counter = 0;
    
    // Save fastq
    char fastq_name[256]; 
    snprintf(fastq_name, sizeof(fastq_name), "%s/%s.fq.gz", output, b.name);
    
    gzFile new_fastq = gzopen(fastq_name, "ab");
    if (!new_fastq) {
        nob_log(NOB_INFO, "Failed to open new fastq_file for appending, exiting");
        return;
    }
    
    // TODO: change this to read_len_max
    char first_part[1000];
    char last_part[1000];
    
    for (size_t j = 0; j < reads.count; ++j) {
        Read r = reads.items[j];
        size_t first_part_start = 0;
        substring(r.seq, first_part_start, barcode_pos, first_part); 
        size_t last_part_start = r.length - barcode_pos;
        substring(r.seq, last_part_start, barcode_pos, last_part); 
        
        // Check for barcode in 5' end
        int match_first_fw = levenshtein_distance(first_part, b.fw, k);
        if (match_first_fw != -1) {
            counter++;
            if (trim) {
                append_read_to_gzip_fastq_trim(new_fastq, &r, match_first_fw, r.length);
            } else {
                append_read_to_gzip_fastq(new_fastq, &r);
            }
            continue;
        }
        
        // Check for reverse barcode in 3' end
        int match_last_rv = levenshtein_distance(last_part, b.fw_comp, k);
        if (match_last_rv != -1) {
            counter++;
            if (trim) {
                append_read_to_gzip_fastq_trim(new_fastq, &r, 0, last_part_start + match_last_rv - b.rv_length);
            } else {
                append_read_to_gzip_fastq(new_fastq, &r);
            }
        }
    }
    gzclose(new_fastq);
    
    pthread_mutex_lock(s_mutex);
    nob_log(NOB_INFO, "barcode: %30s matches: %i", b.name, counter);
    fprintf(S_FILE, "%s,%i\n", b.name, counter);
    pthread_mutex_unlock(s_mutex);
}


void run_nanomux_two(void *arg) {
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

void run_nanomux_one(void *arg) {
    NanomuxData *data = (NanomuxData *)arg;
    Reads local_reads = data->reads;

    process_single_barcode_forward(
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

bool must_be_digit(const char *arg) {
    for(; *arg != '\0'; arg++) {
        if (!isdigit(*arg)) return false;
    }
    return true;
}

char *usage = 
"[USAGE]: nanomux \n"
"   -b <barcode>             Path to barcode file.\n"
"   -f <fastq>               Path to fastq file.\n"   
"   -r <read_len_min>        Minimum length of read.\n"
"   -R <read_len_max>        Maximum length of read.\n"
"   -p <barcode_position>    Position of barcode.\n"
"   -k <mismatches>          Number of misatches allowed.\n"
"   -o <output>              Name of output folder.\n"
"   -t <trim_option>         Trim reads from adapters or not?                      [Options]: trim | notrim.\n"
"   -s <split_option>        Split concatenated reads based on adapter occurance?. [Options]: split | nosplit.\n"
"   -j <threads>             Number of threads to use.                             Default: 1\n";

int main(int argc, char **argv) {    
    
    // CLI arguments
    char *barcode_file = NULL;
    char *fastq_file = NULL;
    int read_len_min = -1;
    int read_len_max = -1;
    size_t barcode_pos = 0;
    int k = -1;
    char *output = NULL;
    char *trim_option = NULL;
    char *split_option = NULL;
    int num_threads = 1;

    bool trim, split;
    
    // parse arguments
    int c;
    while ((c = getopt(argc, argv, "b:f:r:R:p:k:o:t:s:j:")) != -1) {
        switch (c) {
            case 'b':
                barcode_file = optarg;
                break;
            case 'f':
                fastq_file = optarg;
                break;
            case 'r':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-r must be a digit");
                    printf("%s", usage);
                    return 1;
                }
                read_len_min = atoi(optarg);
                break;
            case 'R':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-R must be a digit");
                    printf("%s", usage);
                    return 1;
                }
                read_len_max = atoi(optarg);
                break;
            case 'p':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-p must be a digit");
                    printf("%s", usage);
                    return 1;
                }
                barcode_pos = (size_t)atoi(optarg);
                break;
            case 'k':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-k must be a digit");
                    printf("%s", usage);
                    return 1;
                }
                k = atoi(optarg);
                break;
            case 'o':
                output = optarg;
                break;
            case 't':
                if (strcmp(optarg, "trim") == 0) {
                    trim = true;
                } else if (strcmp(optarg, "notrim") == 0) {
                    trim = false;
                } else {
                    nob_log(NOB_ERROR, "-t must be either 'trim' or 'notrim'");
                    printf("%s", usage);
                    return 1;
                }
                trim_option = optarg;
                break;
            case 's':
                if (strcmp(optarg, "split") == 0) {
                    split = true;
                } else if (strcmp(optarg, "nosplit") == 0) {
                    split = false;
                } else {
                    nob_log(NOB_ERROR, "-s must be either 'split' or 'nosplit'");
                    printf("%s", usage);
                    return 1;
                }
                split_option = optarg;
                break;
            case 'j':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-j must be a digit");
                    printf("%s", usage);
                    return 1;
                }
                num_threads = atoi(optarg);
                break;
            case '?':
                return 1;
        }
    }

    if (!barcode_file || !fastq_file || read_len_min == -1 || read_len_max == -1 || 
        barcode_pos == 0 || k == -1 || !output || !trim_option || !split_option) {
        nob_log(NOB_ERROR, "Error: Missing required arguments\n");
        printf("%s", usage);
        return 1;
    }

    printf("\n");
    nob_log(NOB_INFO, "Running nanomux");
    nob_log(NOB_INFO, "Read len min: %i", read_len_min);
    nob_log(NOB_INFO, "Read len max: %i", read_len_max);
    nob_log(NOB_INFO, "Barcode position: 0 -> %zu", barcode_pos);
    nob_log(NOB_INFO, "k: %i", k);
    nob_log(NOB_INFO, "Trim option: %s", trim_option);
    nob_log(NOB_INFO, "threads: %i", num_threads);
    nob_log(NOB_INFO, "Split option: %s", split_option);
    printf("\n");

    if (!nob_mkdir_if_not_exists(output)) {
        nob_log(NOB_ERROR, "exiting");
        return 1;
    }
   
    // create log file
    char log_file[256];
    snprintf(log_file, sizeof(log_file), "%s/nanomux.log", output);
    FILE *LOG_FILE = fopen(log_file, "ab");
    if (LOG_FILE == NULL) {
        nob_log(NOB_ERROR, "Could create log file");
        return 1;
    }
    
    // create summary file
    char summary_file[256];
    snprintf(summary_file, sizeof(summary_file), "%s/nanomux_matches.csv", output);
    FILE *S_FILE = fopen(summary_file, "ab");
    if (S_FILE == NULL) {
        nob_log(NOB_ERROR, "Could not create summary file");
        fclose(LOG_FILE);
        return 1;
    }
    fprintf(S_FILE, "barcode,matches\n");
    
    fprintf(LOG_FILE, "Nanomux\n\n");
    fprintf(LOG_FILE, "Barcodes: %s\n", barcode_file);
    fprintf(LOG_FILE, "Fastq: %s\n", fastq_file);
    fprintf(LOG_FILE, "Read len min: %i\n", read_len_min);
    fprintf(LOG_FILE, "Read len max: %i\n", read_len_max);
    fprintf(LOG_FILE, "Barcode position: %zu\n", barcode_pos);
    fprintf(LOG_FILE, "K: %i\n", k);
    fprintf(LOG_FILE, "Output folder: %s\n", output);
    fprintf(LOG_FILE, "Trim option: %s\n", trim_option);
    fprintf(LOG_FILE, "Split option: %s\n", split_option);
    
    nob_log(NOB_INFO, "Parsing barcode file %s", barcode_file);
    Nob_String_Builder sb = {0};
    if (!nob_read_entire_file(barcode_file, &sb)) return 1;
    Nob_String_View content = nob_sv_from_parts(sb.items, sb.count);
    Barcodes barcodes = parse_barcodes(content);
    
    // If the number fields in barcode file is 3 --> the file contains fw and rv
    // If the number fields in barcode file is 2 --> the file contains only fw
    int num_fields = num_barcode_fields(barcode_file);
    
    nob_log(NOB_INFO, "Parsing fastq file %s", fastq_file);
    Reads reads = parse_fastq(fastq_file, read_len_min, read_len_max, split);
    nob_log(NOB_INFO, "Number of reads after filtering: %zu", reads.count);
    fprintf(LOG_FILE, "Number of reads after filtering: %zu\n", reads.count);
    printf("\n");

    nob_log(NOB_INFO, "Generating threadpool with %i threads", num_threads);
    printf("\n");
    threadpool thpool = thpool_init(num_threads);
    
    NanomuxDatas nanomux_datas = {0};
    pthread_mutex_t s_mutex = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < barcodes.count; i++) {
        NanomuxData nanomux_data = {
            .barcode = barcodes.items[i],
            .reads = reads,
            .output = output,
            .barcode_pos = barcode_pos,
            .k = k,
            .trim = trim,
            .S_FILE = S_FILE,
            .s_mutex = &s_mutex
        };
        nob_da_append(&nanomux_datas, nanomux_data);
    }
    // forward and reverse barcode
    if (num_fields == 3) {
        nob_log(NOB_INFO, "Barcode file contains forward and reverse barcodes\n");
        for (int i = 0; i < nanomux_datas.count; i++) {
        	thpool_add_work(thpool, run_nanomux_two, (void *)&nanomux_datas.items[i]);
        }
    // only forward barcode
    } else if (num_fields == 2) {
        nob_log(NOB_INFO, "Barcode file only contains one barcode\n");
        for (int i = 0; i < nanomux_datas.count; i++) {
        	thpool_add_work(thpool, run_nanomux_one, (void *)&nanomux_datas.items[i]);
        }
    } else {
        nob_log(NOB_ERROR, "Your barcode file is wrong, please look at the documentation");
        exit(1);
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
    nob_log(NOB_INFO, "nanomux done!");
    
    return 0;
}