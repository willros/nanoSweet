# nanoSweet

## Version
2.0.0

Demultiplex your nanopore (or other) reads! `nanomux` fuzzy matching useful for noisy reads. Written entirely in C.

In the repo, you can also find `nanotrim` – a small threaded program which you can use to filter out reads with a mean quality and between length between min and max. 

Experimental - The repo contains `nanodup` – a small threaded program to deduplicate all the reads and saving information about duplication status.

## Quick Start
```bash
$ git clone https://github.com/willros/nanoSweet.git
$ cd nanoSweet/
$ cc -o nob nob.c
$ ./nob
```

## test nanomux

To get the help message, run `./nanomux`:
```bash
    -b
        Path to barcode file (MANDATORY)
        Default: 
    -f
        Path to fastq file (MANDATORY)
        Default: 
    -o
        Name of output folder (MANDATORY)
        Default: 
    -p
        Position of barcode
        Default: 50
    -k
        Number of mismatches allowed
        Default: 0
    -t
        Trim reads from adapters or not
    -j
        Number of threads to use
        Default: 1
    -help
        Print this help to stdout and exit with 0
    -v
        Print the current version
```

Simple test command:
```bash
./nanomux -b tests/bc_test.csv -f tests/test.fastq -o TEST_NANOMUX -p 100 -k 1 -j 4 -t 
```

The barcode file **MUST** look like this:
```csv
# DUAL BARCODE EXAMPLE
name,forward,reverse
barcode1,ACTATCTACTA,GAGCATGTCGTA
barcode2,AGCGTATGCTGGTA,AGCATGCTATCG

# SINGLE BARCODE EXAMPLE
name,forward
barcode1,ACTATCTACTA
barcode2,AGCGTATGCTGGTA
```

## test nanotrim
To get the help message, run `./nanotrim`:
```bash
    -f
        Path to input folder or file (MANDATORY)
        Default: 
    -o
        Name of output folder (MANDATORY)
        Default: 
    -r
        Minimum read length
        Default: 0
    -R
        Maximum read length
        Default: 1000000
    -q
        Minimum quality
        Default: 0
    -j
        Number of threads to use
        Default: 1
    -help
        Print this help to stdout and exit with 0
    -v
        Print the current version
```


Simple test command:
```bash
./nanotrim -f tests/test.fastq -o TEST_NANOTRIM -r 100 -R 2000 -q 15 
```


## Credit
`nanoSweet` uses `kseq.h` for fastq parsing, and `nob.h`, written by [@tsoding](https://www.github.com/tsoding), for overall useful functions!  
It also uses `thpool.h` by Johan Hanssen Seferidis.

## Change

- 2025-11-07
    - nanomux uses read buffer to process reads now. It does not read all reads into memory anymore.
    - added common.h for shared functions and structures. 




