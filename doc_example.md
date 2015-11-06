---
title: Example &#58; Predicting antibiotic resistance
tags: [getting-started]
keywords: start, introduction, example, kover
last_updated: August 12, 2015
summary: "This page will walk you through an example application of Kover."
---

In this example, we show how Kover can be applied to genomic data in order to obtain an interpretable model of a phenotype.

## Example data

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

## Splitting the dataset

In order to measure the accuracy of the model obtained using Kover, we must split the dataset into a training set and a 
testing set. The training set will be used to learn a model and the testing set will be used to estimate its accuracy.
A Kover dataset can contain multiple different splits of the data.

Moreover, Kover being a machine learning algorithm, it has hyperparameters, i.e. user specified parameters, that must
be adjusted based on the data. Kover can use [multiple strategies](todo.com) for choosing the hyperparameter values. In this example, we
will use [k-fold cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)#k-fold_cross-validation), 
which is the most commonly used hyperparameter selection strategy in machine learning.

The following command creates a split of the data called "example_split", which uses 2/3 of the data for training and
1/3 for testing. It also creates 5 cross-validation folds. The data splitting is done randomly with random seed 42.

```
kover dataset split --dataset example.kover --id example_split --train-size 0.666 --folds 5 --random-seed 42 --progress
```

## Learning a model

Now that we have created and split the dataset, we are ready to learn a Rifampicin resistance model.