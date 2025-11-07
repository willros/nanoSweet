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
Run `./nanomux` to get the help message.

## nanotrim
Run `./nanotrim` to get the help message.


## Credit
`nanoSweet` uses `kseq.h` for fastq parsing, and `nob.h`, written by [@tsoding](https://www.github.com/tsoding), for overall useful functions!  
It also uses `thpool.h` by Johan Hanssen Seferidis.

## Change

- 2025-11-07
    - nanomux uses read buffer to process reads now. It does not read all reads into memory anymore.
    - added common.h for shared functions and structures. 




