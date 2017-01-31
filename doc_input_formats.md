---
title: Input Data Formats
keywords: kover, machine learning, data, format, input, reads, contigs, k-mer
last_updated: January 31, 2017
tags:
summary: "An overview of the file formats accepted by Kover."
---

## Genomic Data

Kover currently accepts genomic data in three formats:

* reads: a set of FASTQ files containing genomic reads

* contigs: a set of FASTA files containing assembled genomic sequences

* k-mer matrix: a matrix giving the presence/absence of each k-mer in each genome


### Reads

In this case, the genomic data is available as a set of [FASTQ files](https://en.wikipedia.org/wiki/FASTQ_format).
There can be more than one read file per genome.
You must provide a tab-separated value (TSV) file relating each genome to a folder containing its reads files.
It should have the following format:

| | |
| :-: | :-: | :-: | :-: | :-: |
|GenomeID_1| Read_folder_1|
|GenomeID_2| Read_folder_1|
| ...      | ... |
|GenomeID_m| Read_folder_m|

**Please make sure that the genome identifiers in the TSV file match the ones in the metadata.*


### Contigs

In this case, the genomic data is available as a set of [FASTA files](https://en.wikipedia.org/wiki/FASTA_format) (one per genome).
Each file contains a set of contigs, which are assembled genomic sequences. 
You must provide a tab-separated value (TSV) file relating each FASTA file to a genome.
It should have the following format:

| | |
| :-: | :-: | :-: | :-: | :-: |
|GenomeID_1| FASTA_Path_1|
|GenomeID_2| FASTA_Path_2|
| ...      | ... |
|GenomeID_m| FASTA_Path_m|

**Please make sure that the genome identifiers in the TSV file match the ones in the metadata.*


### K-mer matrix

In this case, the genomic data is available as a tab-separated value (TSV) file where:

* The first line is a header, with the first column labelled "kmers" and the remaining columns labelled with
genome identifiers. For example, for a study based on 100 genomes, there should be 101 columns in the file.

* Each of the remaining lines gives the presence or absence of a k-mer in each genome. Each line starts with the
k-mer sequence and the remaining columns contain a 0 if the k-mer is absent in the genome or a 1 if it is
present.

| kmers | GenomeID_1 | GenomeID_2 | ... | GenomeID_m |
| :-: | :-: | :-: | :-: | :-: |
| kmer_1 | 1 | 0 | ... | 0 |
| kmer_2 | 0 | 1 | ... | 1 |
| ... | ... | ... | ... | ... |
| kmer_n | 1 | 0 | ... | 1 |

Such a matrix can be generated with [Ray Surveyor](https://github.com/zorino/RaySurveyor-Tutorial).

**Please make sure that the genome identifiers in the k-mer matrix match the ones in the metadata.*


## Metadata

The metadata must be provided as a two-column TSV file. Each line contains a genome identifier and a binary value (0 or
1) indicating its associated phenotype. The meaning of each value is arbitrary, as it only specifies a grouping of
genomes. This is the phenotypic data that will be used to train the learning algorithm.

| | |
| --- | :-: |
| GenomeID_1 | 1 |
| GenomeID_2 | 0 |
| ... | ... |
| GenomeID_m | 0 |

**Notice that there is no header.*

**Please make sure that the genome identifiers in the metadata match the ones in the genomic data.*