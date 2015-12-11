---
title: Kover dataset tools
tags: 
keywords: dataset, tools, kover, genomics, k-mer, machine learning
last_updated: December 11, 2015
summary: "Overview of the Kover dataset tools"
---

## What is a Kover dataset?

Coming soon...

## Creating datasets

```
usage: kover dataset create [-h] --genome-type {tsv} --genome-source
                            GENOME_SOURCE --output OUTPUT
                            [--phenotype-name PHENOTYPE_NAME]
                            [--phenotype-metadata PHENOTYPE_METADATA]
                            [--compression COMPRESSION] [-x] [-v]

Creates a Kover dataset from genomic data and optionally phenotypic metadata

optional arguments:
  -h, --help            show this help message and exit
  --genome-type {tsv}   The type of source for the genomic data.
  --genome-source GENOME_SOURCE
                        The source of the genomic data.
  --output OUTPUT       The output Kover dataset.
  --phenotype-name PHENOTYPE_NAME
                        The name of the phenotype.
  --phenotype-metadata PHENOTYPE_METADATA
                        A file containing the metadata.
  --compression COMPRESSION
                        The gzip compression level (0 - 9). 0 means no
                        compression. The default value is 4.
  -x, --progress        Shows a progress bar for the execution.
  -v, --verbose         Sets the verbosity level.
```

## Splitting datasets

```
usage: kover dataset split [-h] --dataset DATASET --id ID
                           [--train-size TRAIN_SIZE]
                           [--random-seed RANDOM_SEED] [--folds FOLDS] [-v]
                           [-x]

Splits a kover dataset file into a training set, a testing set and optionally
cross-validation folds

optional arguments:
  -h, --help            show this help message and exit
  --dataset DATASET     The Kover dataset to split.
  --id ID               The identifier of the split.
  --train-size TRAIN_SIZE
                        The proportion of the data used for training (default
                        is 0.5).
  --random-seed RANDOM_SEED
                        The random seed used for the split (If not provided a
                        random value between 0 and 4294967295 will be used).
  --folds FOLDS         The number of cross-validation folds (default is 0 for
                        none, the minimum value is 2).
  -v, --verbose         Sets the verbosity level.
  -x, --progress        Shows a progress bar for the execution.
  ```
