---
title: Data Manipulation Module
tags: 
keywords: dataset, tools, kover, genomics, k-mer, machine learning
last_updated: January 31, 2017
summary: "Overview of the data manipulation utility"
---

## Creating a dataset

This command is used to combine genomic and phenotypic data into a Kover dataset.

### From reads

Use this command if the genomic data consists of unassembled genomic reads (see [input data format](doc_input_formats.html)).
The [DSK k-mer counter](https://gatb.inria.fr/software/dsk/) is used to count the k-mers present in the genomes.

```
usage: kover dataset create from-reads [-h] --genomic-data GENOMIC_DATA
                                       [--phenotype-description PHENOTYPE_DESCRIPTION]
                                       [--phenotype-metadata PHENOTYPE_METADATA]
                                       --output OUTPUT [--kmer-size KMER_SIZE]
                                       [--kmer-min-abundance KMER_MIN_ABUNDANCE]
                                       [--singleton-kmers] [--n-cpu N_CPU]
                                       [--compression COMPRESSION]
                                       [--temp-dir TEMP_DIR] [-x] [-v]

Creates a Kover dataset from genomic data and optionally phenotypic metadata

optional arguments:
  -h, --help            show this help message and exit
  --genomic-data GENOMIC_DATA
                        A tab-separated file with one line per genome in the
                        format GENOME_ID{tab}PATH, where the path refers to a
                        directory containing the genome's reads in fastq(.gz)
                        files.
  --phenotype-description PHENOTYPE_DESCRIPTION
                        An informative description that is assigned to the
                        phenotypic metadata.
  --phenotype-metadata PHENOTYPE_METADATA
                        A file containing the phenotypic metadata.
  --output OUTPUT       The Kover dataset to be created.
  --kmer-size KMER_SIZE
                        The k-mer size (max is 128). The default is 31.
  --kmer-min-abundance KMER_MIN_ABUNDANCE
                        The minimum number of times a k-mer must be found in a
                        read file in order to be considered. All k-mers that
                        do not meet this threshold are discarded. This value
                        should be chosen based on genome coverage (ex: 100x
                        coverage -> you could use 10). The default is 1.
  --singleton-kmers     Include k-mers that only occur in one genome. Disabled
                        by default.
  --n-cpu N_CPU, --n-cores N_CPU
                        The number of cores used by DSK. The default value is
                        0 (all cores).
  --compression COMPRESSION
                        The gzip compression level (0 - 9). 0 means no
                        compression. The default value is 4.
  --temp-dir TEMP_DIR   Output directory for temporary files. The default is
                        the system's temp dir.
  -x, --progress        Shows a progress bar for the execution.
  -v, --verbose         Sets the verbosity level.
```

### From contigs

Use this command if the genomic data consists of assembled genomes (one fasta file per genome) (see [input data format](doc_input_formats.html)).
The [DSK k-mer counter](https://gatb.inria.fr/software/dsk/) is used to count the k-mers present in the genomes.

```
usage: kover dataset create from-contigs [-h] --genomic-data GENOMIC_DATA
                                         [--phenotype-description PHENOTYPE_DESCRIPTION]
                                         [--phenotype-metadata PHENOTYPE_METADATA]
                                         --output OUTPUT
                                         [--kmer-size KMER_SIZE]
                                         [--singleton-kmers] [--n-cpu N_CPU]
                                         [--compression COMPRESSION]
                                         [--temp-dir TEMP_DIR] [-x] [-v]

Creates a Kover dataset from genomic data and optionally phenotypic metadata

optional arguments:
  -h, --help            show this help message and exit
  --genomic-data GENOMIC_DATA
                        A tab-separated file with one line per genome in the
                        format GENOME_ID{tab}PATH, where the path refers to a
                        fasta file containing the genome's contigs.
  --phenotype-description PHENOTYPE_DESCRIPTION
                        An informative description that is assigned to the
                        phenotypic metadata.
  --phenotype-metadata PHENOTYPE_METADATA
                        A file containing the phenotypic metadata.
  --output OUTPUT       The Kover dataset to be created.
  --kmer-size KMER_SIZE
                        The k-mer size (max is 128). The default is 31.
  --singleton-kmers     Include k-mers that only occur in one genome. Disabled
                        by default.
  --n-cpu N_CPU, --n-cores N_CPU
                        The number of cores used by DSK. The default value is
                        0 (all cores).
  --compression COMPRESSION
                        The gzip compression level (0 - 9). 0 means no
                        compression. The default value is 4.
  --temp-dir TEMP_DIR   Output directory for temporary files. The default is
                        the system's temp dir.
  -x, --progress        Shows a progress bar for the execution.
  -v, --verbose         Sets the verbosity level.
```

### From a k-mer matrix
Use this command if the genomic data consists of a matrix giving the presence or absence of each k-mer in each genome (see [input data format](doc_input_formats.html)).

```
usage: kover dataset create from-tsv [-h] --genomic-data GENOMIC_DATA
                                     [--phenotype-description PHENOTYPE_DESCRIPTION]
                                     [--phenotype-metadata PHENOTYPE_METADATA]
                                     --output OUTPUT
                                     [--compression COMPRESSION] [-x] [-v]

Creates a Kover dataset from genomic data and optionally phenotypic metadata

optional arguments:
  -h, --help            show this help message and exit
  --genomic-data GENOMIC_DATA
                        A tab-separated file containing the k-mer matrix.
  --phenotype-description PHENOTYPE_DESCRIPTION
                        An informative description that is assigned to the
                        phenotypic metadata.
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
                           [--train-size TRAIN_SIZE] [--train-ids TRAIN_IDS]
                           [--test-ids TEST_IDS] [--folds FOLDS]
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
                        Alternatively, you can specify which genomes to use
                        for training and testing by using --train-ids and
                        --test-ids.
  --train-ids TRAIN_IDS
                        File containing the identifiers of the genomes used to
                        train the learning algorithm. If you provide a value
                        for this argument, you must also provide a value for
                        --test-ids. File format: one id per line
  --test-ids TEST_IDS   File containing the identifiers of the genomes used to
                        evaluate the accuracy of the model generated. If you
                        provide a value for this argument, you must also
                        provide a value for --train-ids. File format: one id
                        per line
  --folds FOLDS         The number of k-fold cross-validation folds to create
                        (default is 0 for none, the minimum value is 2). Folds
                        are required for using k-fold cross-validation in
                        'kover learn'.
  --random-seed RANDOM_SEED
                        A random seed used for randomly splitting the data. A
                        specific seed will always lead to the same split. If
                        not provided, it is set randomly.
  -v, --verbose         Sets the verbosity level.
  -x, --progress        Shows a progress bar for the execution
```

## Listing information about a dataset

This command is used to list any information about a Kover dataset.

```
usage: kover dataset info [-h] --dataset DATASET [--all] [--genome-type]
                          [--genome-source] [--genome-ids] [--genome-count]
                          [--kmers] [--kmer-len] [--kmer-count]
                          [--phenotype-description] [--phenotype-metadata]
                          [--phenotype-tags] [--splits] [--uuid]
                          [--compression] [--classification-type]

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
  --phenotype-description
                        Prints the description that was assigned to the
                        phenotype.
  --phenotype-metadata  Prints the path of the file from which the phenotypic
                        metadata was acquired.
  --phenotype-tags      Prints the phenotype tags associated to the dataset
  --splits              Prints the lists of splits of the dataset that are
                        available for learning.
  --uuid                Prints the unique identifier of the dataset.
  --compression         Prints the data compression options of the dataset.
  --classification-type
                        Prints the dataset classification type.
```
