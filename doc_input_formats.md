---
title: Input Data Formats
keywords: kover, machine learning, data, format, input
last_updated: April 14, 2016
tags:
summary: "An overview of the file formats accepted by Kover."
---

## Genomic Data

Kover currently accepts genomic data in a specific tab-separated value (TSV) file format, which corresponds to the k-mer
matrix produced by [Ray Surveyor](https://github.com/zorino/RaySurveyor-Tutorial).

* The first line should be a header, with the first column labelled "kmers" and the remaining columns labelled with
genome identifiers. For example, for a study based on 100 genomes, there should be 101 columns in the file.

* Each of the remaining lines gives the presence or absence of a k-mer in each genome. Each line should start with the
k-mer sequence and the remaining columns should contain a 0 if the k-mer is absent in the genome or a 1 if it is
present.

| kmers | GenomeID_1 | GenomeID_2 | ... | GenomeID_m |
| :-: | :-: | :-: | :-: | :-: |
| kmer_1 | 1 | 0 | ... | 0 |
| kmer_2 | 0 | 1 | ... | 1 |
| ... | ... | ... | ... | ... |
| kmer_n | 1 | 0 | ... | 1 |

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

**Please make sure that the genome identifiers in the metadata match the ones in the k-mer matrix.*