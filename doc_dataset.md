---
title: Data Manipulation Module
tags: 
keywords: dataset, tools, kover, genomics, k-mer, machine learning
last_updated: March 27, 2016
summary: "Overview of the data manipulation utilities"
---

## Creating a dataset

This command is used to combine genomic and phenotypic data into a Kover dataset.

```
usage: kover dataset create [-h] --genome-type {tsv} --genome-source
                            GENOME_SOURCE [--phenotype-name PHENOTYPE_NAME]
                            [--phenotype-metadata PHENOTYPE_METADATA] --output
                            OUTPUT [--compression COMPRESSION] [-x] [-v]

Creates a Kover dataset from genomic data and optionally phenotypic metadata

optional arguments:
  -h, --help            show this help message and exit
  --genome-type {tsv}   The format in which the genomic data is provided. See
                        documentation for details.
  --genome-source GENOME_SOURCE
                        The genomic data.
  --phenotype-name PHENOTYPE_NAME
                        An informative name that is assigned to the phenotypic
                        metadata.
  --phenotype-metadata PHENOTYPE_METADATA
                        A file containing the phenotypic metadata.
  --output OUTPUT       The Kover dataset to be created.
  --compression COMPRESSION
                        The gzip compression level (0 - 9). 0 means no
                        compression. The default value is 4.
  -x, --progress        Shows a progress bar for the execution.
  -v, --verbose         Sets the verbosity level.
```

## Splitting a dataset

This command is used to split a Kover dataset into a training set, a testing set and optionally cross-validation folds.
This must be done prior to learning models from the data.

```
usage: kover dataset split [-h] --dataset DATASET --id ID
                           [--train-size TRAIN_SIZE] [--folds FOLDS]
                           [--random-seed RANDOM_SEED] [-v] [-x]

Splits a kover dataset file into a training set, a testing set and optionally
cross-validation folds

optional arguments:
  -h, --help            show this help message and exit
  --dataset DATASET     The Kover dataset to be split.
  --id ID               A unique identifier that will be assigned to the
                        split.
  --train-size TRAIN_SIZE
                        The proportion of the data that will be reserved for
                        training the learning algorithm (default is 0.5).
  --folds FOLDS         The number of k-fold cross-validation folds to create
                        (default is 0 for none, the minimum value is 2). Folds
                        are required for using k-fold cross-validation in
                        'kover learn'.
  --random-seed RANDOM_SEED
                        A random seed used for randomly splitting the data. A
                        specific seed will always lead to the same split. If
                        not provided, it is set randomly.
  -v, --verbose         Sets the verbosity level.
  -x, --progress        Shows a progress bar for the execution.
```

## Listing information about a dataset

This command is used to list any information about a Kover dataset.

```
usage: kover dataset info [-h] --dataset DATASET [--all] [--genome-type]
                          [--genome-source] [--genome-ids] [--genome-count]
                          [--kmers] [--kmer-len] [--kmer-count]
                          [--phenotype-name] [--phenotype-metadata] [--splits]
                          [--uuid] [--compression]

Prints information about the content of a dataset

optional arguments:
  -h, --help            show this help message and exit
  --dataset DATASET     The Kover dataset for which you require information.
  --all                 Prints all the available information.
  --genome-type         Prints the type of genomic data that was used to
                        create the dataset.
  --genome-source       Prints the source (e.g.: path) from which the genomic
                        data was acquired.
  --genome-ids          Prints the identifiers of the genomes in the dataset.
  --genome-count        Prints the number of genomes in the dataset.
  --kmers               Prints the sequence of each k-mer in the dataset
                        (fasta).
  --kmer-len            Prints the length of the k-mers in the dataset.
  --kmer-count          Prints the number of k-mers in the dataset.
  --phenotype-name      Prints the identifier that was assigned to the
                        phenotype.
  --phenotype-metadata  Prints the path of the file from which the phenotypic
                        metadata was acquired.
  --splits              Prints the lists of splits of the dataset that are
                        available for learning.
  --uuid                Prints the unique identifier of the dataset.
  --compression         Print the data compression options of the dataset.
```