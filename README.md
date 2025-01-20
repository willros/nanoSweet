# nanoSweet

Demultiplex your nanopore (or other) reads! `nanomux` fuzzy matching useful for noisy reads. Written entirely in C, so should compile on most systems. The repo contains a windows branch that uses something else than `pthreads` for threading.

In the repo, you can also find `nanotrim` – a small threaded program which you can use to filter out reads with a mean quality and between length between min and max. 

The repo also contains `nanodup` – a small threaded program to deduplicate all the reads and saving information about duplication status.

## Quick Start
```bash
$ git clone https://github.com/willros/nanomux_c.git
$ cd nanomux_c/
$ cc -o nob nob.c
$ ./nob
```

## test nanomux
```bash
$ ./nanomux -b tests/bc_test.csv -f tests/test.fastq -r 600 -R 2000 -p 200 -k 1 -o new_nanomux -t trim -s split -j 4
```

## test nanotrim
```bash
$ ./nanotrim -i tests/test.fastq -r 2000 -R 10000 -q 20 -t 4 -o test_nanotrim
```

## test nanodup
```bash
$ ./nanodup -i tests/test.fastq -o test_nanodup
```

## nanomux
Run `./nanomux` to get the help message:
```bash
[USAGE]: nanomux 
   -b <barcode>             Path to barcode file.
   -f <fastq>               Path to fastq file.
   -r <read_len_min>        Minimum length of read.
   -R <read_len_max>        Maximum length of read.
   -p <barcode_position>    Position of barcode.
   -k <mismatches>          Number of misatches allowed.
   -o <output>              Name of output folder.
   -t <trim_option>         Trim reads from adapters or not?                      [Options]: trim | notrim.
   -s <split_option>        Split concatenated reads based on adapter occurance?. [Options]: split | nosplit.
   -j <threads>             Number of threads to use.                             Default: 1
```
* `barcode_pos`: where to search for the barcode in the ends. if barcode_pos == 200 and the read length is 1000, the barcodes will be searched for from position 0 -> 200 and 800 -> 1000.
* `k`: allowed number of mismatches
* `trim_option`: To trim the reads to the left and right of the found barcode.
* `split_option`: Try to split reads longer than 2500 if an adapter sequence is found in the middle. One splitted read results in two new reads, with the suffix `_1` and `_2`.

`barcode_file.csv` **MUST** have the follwing shape:
* columns:
    * name
    * forward barcode
    * reverse barcode

It must look like the example below. 
```csv
bc1,ATACGATGCTA,GTCGATGTCTGA
bc2,GACACACAC,GTCGATTGATG
```

## nanotrim
Run `./nanotrim` to get the help message:
```bash
[USAGE]: nanotrim -i <input> [options]
"[USAGE]: nanotrim -i <input> [options]\n"
"   -i    <input>             Path of folder or file\n"
"   -o    <output>            Name of output folder.\n"
"   -r    <read_length_min>   Minium length of read.    Optional: Default 1\n"
"   -R    <read_length_max>   Minium length of read.    Optional: Default INT_MAX\n"
"   -q    <quality>           Minimum quality of read.  Optional: Default 1\n"
"   -t    <threads>           Number of threads to use. Optional: Default 1\n";
```

`nanotrim` saves the trimmed reads to a fastq in the specified output folder. Uses same name as the original file, but with the suffix `.filtered`. The input can be a single file or an entire folder – `nanotrim` knows can distinguish between the two.  


Example of output:
```bash
$ ./nanotrim -i tests/test.fastq -r 2000 -R 10000 -q 20 -t 4
$ tests/test.fastq: 4788 raw reads (29 passed) --> To short: 4169  | To long: 3     | Low quality: 587
```

`nanotrim` also produces a log file with the above information in a `.csv` format, saved in the output path.

```csv
file,raw_reads,passed_reads,short,long,bad_quality
tests/test.fastq,4788,29,4169,3,587
```

## nanodup
Run `./nanodup` to get the help message:
```bash
[ERROR] You must provide the input path
[ERROR] [USAGE]: nanodup -i <input> -o <output> [options]
   -i    <input>             Path of folder or file
   -o    <output>            Name of output folder.
   -t    <threads>           Number of threads to use. Optional: Default 1
```

`nanodup` keeps removes all the duplicated reads, *i.e.* identical reads. It only saves the first read if the read is duplicated. `nanodup` also produces a log file with the information about the reads and saves all the de-duplicated reads to a new fastq file. 

Example of output:
```bash
$ ./nanodup -i tests/test.fastq -o test_nanodup
$ [INFO] tests/test.fastq contained: 0 duplicates
```


## Credit
`nanomux_c` uses `kseq.h` for fastq parsing, and `nob.h`, written by [@tsoding](https://www.github.com/tsoding), for overall useful functions!  
It also uses `thpool.h` by Johan Hanssen Seferidis.

## TODO
- [x] trim barcodes
- [x] multi threading
- [x] read splitting
- [x] Change to threadpool in nanomux
- [x] Add logging to nanotrim
- [] Fix the windows branch with gzip append and nanodup




