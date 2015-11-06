---
title: Example &#58; Predicting antibiotic resistance
tags: [getting-started]
keywords: start, introduction, example, kover
last_updated: August 12, 2015
summary: "This page will walk you through an example application of Kover."
---

In this example, we show how Kover can be applied to genomic data in order to obtain an interpretable model of a phenotype.

## Data

You will need to download the [example data](blob) (~180 Mb). The data contains the genomes of 141 *Mycobacterium 
tuberculosis* isolates, along with their susceptibility to Rifampicin, an antibiotic.

The genomes were assembled using the [SPades](http://bioinf.spbau.ru/spades) genome assembler and subsequently split into
k-mers of length 31 using [Ray Surveyor](https://github.com/zorino/RaySurveyor-Tutorial). The tsv-format matrix was produced
by Ray Surveyor.

The genomes and the metadata were obtained from: Merker, Matthias, et al. "Evolutionary history and global spread of the Mycobacterium tuberculosis Beijing lineage." *Nature genetics* 47.3 (2015): 242-249.

## Creating a dataset

Before learning a Rifampicin resistance model from the example data, we must package the genomic and phenotypic data into a [Kover Dataset](). 
To convert the example data into such a dataset, use the following command:

```
kover dataset create --genome-type tsv --genome-source KmerMatrix.tsv --phenotype-name "Rifampicin resistance" --phenotype-metadata metadata.tsv --output example.kover --progress
```

This produces a Kover dataset called "example.kover". From now on, you no longer need the original example data files.

You can now use the [kover dataset info]() command to print information about the dataset. For example, to list the identifiers
of the genomes contained in the dataset, use:

```
kover dataset info --dataset example.kover --genome-ids
```

## Learning a model