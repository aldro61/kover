---
title: Example &#58; Predicting antibiotic resistance
tags: [getting-started]
keywords: start, introduction, example, kover
last_updated: November 9, 2015
summary: "This page will walk you through an example application of Kover."
---

In this example, we show how Kover can be applied to genomic data in order to obtain an interpretable model of a phenotype.

## Example data

You will need to download the [example data](http://graal.ift.ulaval.ca/adrouin/kover-example-data.zip) (~180 Mb). The data contains the genomes of 141 *Mycobacterium 
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

You can now use the [kover dataset info](todo) command to print information about the dataset. For example, to list the identifiers
of the genomes contained in the dataset, use:

```
kover dataset info --dataset example.kover --genome-ids
```

and to print the number of k-mers in the dataset, use:

```
kover dataset info --dataset example.kover --kmer-count
```

Your dataset contains 9 701 935 k-mers!

## Splitting the dataset

In order to measure the accuracy of the model obtained using Kover, we must split the dataset into a training set and a 
testing set. The training set will be used to learn a model and the testing set will be used to estimate its accuracy.
A Kover dataset can contain multiple splits of the data. The command for splitting a dataset is [kover dataset split](todo).

Moreover, Kover being a machine learning algorithm, it has [hyperparameters](todo), i.e. user specified parameters, that must
be adjusted based on the data. Kover can use [multiple strategies](todo) for choosing the hyperparameter values. In this example, we
will use [k-fold cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)#k-fold_cross-validation), 
which is the most commonly used hyperparameter selection strategy in machine learning.

The following command creates a split of the data called "example_split", which uses 2/3 of the genomes for training and
1/3 for testing. It also creates 5 cross-validation folds. The data splitting is done randomly with 42 as the random seed.

```
kover dataset split --dataset example.kover --id example_split --train-size 0.666 --folds 5 --random-seed 62 --progress
```

## Learning a model

Now that we have created and splitted the dataset, we are ready to learn a Rifampicin resistance model. The [kover learn]()
command is used for learning. The following command tells Kover to learn a model containing at most 5 rules, to try both
conjunction and disjunction models (see [model types](todo)) and the values 0.1, 1.0 and 10.0 for the *p* 
hyperparameter (see [hyperparameters](todo)), while using cross-validation as the hyperparameter selection strategy. 
Moreover, it distributes the cross-validation on 2 CPUs.

```
kover learn --dataset example.kover --split example_split --model-type conjunction disjunction --p 0.1 1.0 10.0 --max-rules 5 --hp-choice cv --n-cpu 2 --progress
```

Kover then uses the obtained model to predict the phenotype of the genomes in the testing set and computes various metrics.
For this example, the obtained model is:

```
Absence(CCCAGCGCCGACAGTCGGCGCTTGTGGGTCA) [Importance: 0.97]
OR
Absence(CGCAACAAGTCAGCGTCCCTGAGGGGGGGCA) [Importance: 0.36]
```

meaning that if any of these sequences is not present in the genome, then the isolate can be considered resistant to Rifampicin.
Notice the simplicity and interpretability of the obtained model. 
The testing set metrics for this model are:

```
Error Rate: 0.04255
Sensitivity: 0.95652
Specificity: 0.95833
Precision: 0.95652
Recall: 0.95652
F1 Score: 0.95652
True Positives: 22.0
True Negatives: 23.0
False Positives: 1.0
False Negatives: 1.0
```

## Subsequent analysis of the model

You could use [Nucleotide Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch) to identify the genomic regions corresponding the sequences targeted by the obtained model. 
In this case, the sequence [GCGCCGACAGTCGGCGCTTGTGGGTCAACCC](https://www.ncbi.nlm.nih.gov/nucleotide/746590776?from=76&to=106) corresponds to the *rpoB* gene, which encodes the RNA polymerase
beta subunit. This k-mer falls exactly in the Rifampicin resistance determining region of the gene. Moreover, the fact
that the selected rule is an absence rule suggests that there are many variant sequences at this position that confer
resistance to Rifampicin. An absence rule can concisely regroup many presence rules.

## Predicting with the obtained model

You could use any k-mer counting tool, such as [Jellyfish](https://github.com/gmarcais/Jellyfish), to extract the k-mers present in a FASTA file and then, easily apply the model.
This feature is currently not included in Kover, but will be added in future releases.
