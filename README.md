# California leaf-nosed bat (*Macrotus californicus*) genome assembly

# Getting started

## Dependencies

 - FastQC (0.11.2 or greater)
 - Snakemake (depends on Python 3)
 - Make

## Configuration

There is a `conf.py` file containing constants and configuration.
Set it up accordingly for your system.

## Running

The Makefile contains the 'bootstrap' for Snakemake (download raw data and prepare
the raw data file list).

```bash
$ make
```

For subsequent runs and other steps you can invoke Snakemake directly:

```bash
$ snakemake <rule name or target file>
```

# Analysis steps

## Raw reads quality control

FastQC generates a report for both Illumina and PacBio raw reads.
The following rule executes FastQC with the appropriate files and publish the reports
to a remote SSH host.

```bash
$ snakemake publish_fastqc
```
