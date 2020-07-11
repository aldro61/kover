---
title: An introduction to Kover with Set Covering Machines
tags: [getting-started]
keywords: learning, set covering machine, scm, model
last_updated: August 8, 2017
summary: "This tutorial will show how to use Kover with the Set Covering Machine algorithm."
---

In this example, we show how the Set Covering Machine algorithm, implemeted in Kover, can be applied to genomic data in order to obtain an interpretable model of a phenotype.
Specifically, we will learn a model that predicts rifampicin resistance in *Mycobacterium tuberculosis*.

## Example data

First, download the [example data](https://github.com/aldro61/kover-tutorial-data/releases/download/1.0.0/kover-example-data.zip) (~180 Mb), which contains the genome of 141 *Mycobacterium
tuberculosis* isolates, along with their susceptibility to rifampicin. There should be two files: KmerMatrix.tsv and metadata.tsv.

The genomes were assembled using the [SPades](http://bioinf.spbau.ru/spades) and split into k-mers of length 31 using
[Ray Surveyor](https://github.com/zorino/RaySurveyor-Tutorial).

The raw data were obtained from: Merker, Matthias, et al. "Evolutionary history and global spread of the Mycobacterium tuberculosis Beijing lineage." *Nature genetics* 47.3 (2015): 242-249.

![#1589F0](https://placehold.it/10/1589F0/000000?text=+) **Note:** Kover also work on reads and contigs (BAM, FASTA, FASTQ, etc.) ([see here for details](doc_input_formats.html)).

![#1589F0](https://placehold.it/10/1589F0/000000?text=+) **Additional links:** [Ray Surveyor Tutorial](https://github.com/zorino/RaySurveyor-Tutorial)


## Creating a dataset

Before using Kover to learn a model, we must package the genomic and phenotypic data into a [Kover dataset](doc_dataset.html#creating-a-dataset), which relies on the HDF5 library to store a compressed representation of the data ([details here](https://github.com/aldro61/kover/wiki/Kover-Dataset-Format)).

To create a dataset, use the following command:

```
kover dataset create from-tsv --genomic-data KmerMatrix.tsv --phenotype-description "Rifampicin resistance" --phenotype-metadata metadata.tsv --output example.kover --progress
```

This produces a dataset file called "example.kover". From now on, you no longer need the original data files.


## Exploring the dataset

You can use the [kover dataset info](doc_dataset.html#listing-information-about-a-dataset) command to print information about the dataset. For example, to list the identifiers
of the genomes contained in the dataset, use:

```
kover dataset info --dataset example.kover --genome-ids
```

To print the number of genomes and k-mers in the dataset, use:

```
kover dataset info --dataset example.kover --genome-count --kmer-count
```

Your dataset contains **141 genomes vs 9 701 935 k-mers**! This is know as the *fat data* setting, which is very different from the *big data* setting in which the number of examples (genomes) is greater than the number of features (k-mers).


## Splitting the dataset

In order to measure the accuracy of the model obtained using Kover, we must split the dataset into a training set and a 
testing set. The training set will be used to learn a model and the testing set will be used to estimate its accuracy.
A Kover dataset can contain multiple splits of the data. The command used for splitting a dataset is [kover dataset split](doc_dataset.html#splitting-a-dataset).

The Set Covering Machine algorithm has [hyperparameters](doc_learning.html#understanding-the-hyperparameters),
which are free parameters that control its behavior and that must be tuned to the data. The most widely used method for setting hyperparameter values
is [k-fold cross-validation](doc_learning.html#k-fold-cross-validation).
In this example, we will use 5-fold cross-validation.

The following command creates a split of the data called "example_split", which uses 2/3 of the genomes for training and
1/3 for testing. It also creates 5 cross-validation folds. The data are partitioned randomly, using 72 as the random seed.

```
kover dataset split --dataset example.kover --id example_split --train-size 0.666 --folds 5 --random-seed 72 --progress
```

## Learning a model

Now that we have created and split the dataset, we are ready to learn a predictive model of Rifampicin resistance in *Mycobacterium tuberculosis*.
The [kover learn](doc_learning.html#command-line-interface) command is used to learn models.
The following command tells Kover to use the Set Covering Machine algorithm to learn a model containing at most 5 rules, to try both
conjunction (logical-AND) and disjunction (logical-OR) models and the values 0.1, 1.0 and 10.0 for the *p*
hyperparameter (see [hyperparameters](doc_learning.html#understanding-the-hyperparameters)), while using cross-validation as the [hyperparameter selection strategy](doc_learning.html#hyperparameter-selection-strategies).
Moreover, it distributes the cross-validation on 2 CPUs.

```
kover learn scm --dataset example.kover --split example_split --model-type conjunction disjunction --p 0.1 1.0 10.0 --max-rules 5 --hp-choice cv --n-cpu 2 --progress
```

Kover then uses the obtained model to predict the phenotype of the genomes in the testing set and computes various metrics.
For this example, the obtained model is:

```
Absence(GCGCCGACAGTCGGCGCTTGTGGGTCAACCC) [Importance: 0.93, 6 equivalent rules]
OR
Absence(ACCAGAACAACCCGCTGTCGGGGTTGACCCA) [Importance: 0.09, 1 equivalent rules]
```

The model indicates that, if any of these sequences is not present in the genome, then the isolate is resistant to rifampicin.
Notice the simplicity and interpretability of the obtained model. 

The testing set metrics for this model are:

```
Error Rate: 0.0
Sensitivity: 1.0
Specificity: 1.0
Precision: 1.0
Recall: 1.0
F1 Score: 1.0
True Positives: 15.0
True Negatives: 32.0
False Positives: 0.0
False Negatives: 0.0
```

![#1589F0](https://placehold.it/10/1589F0/000000?text=+) **Note:**  The randomness used to split the dataset can vary based on the computer and operating system. This could explain
why this example gives slightly different results on your computer. The number of rules in the model, the k-mer sequences
and the accuracy could be slightly different.


## Subsequent analysis of the model

We can now further analyse our model and try to elucidate the nature of the k-mers used by the model. To achieve this, we can use [Nucleotide BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch).

![#c5f015](https://placehold.it/10/c5f015/000000?text=+) **Example:** Use BLAST to search for the following k-mer, which is of high importance in our model:

```
GCGCCGACAGTCGGCGCTTGTGGGTCAACCC
```

You should find that the sequence maps to the *rpoB* gene, which encodes the RNA polymerase
beta subunit ([see here](https://www.ncbi.nlm.nih.gov/nucleotide/746590776?from=76&to=106)). Interestingly, this k-mer is in the rifampicin resistance determining region of the gene, so it seems like we have successfully identified a resistance determinant, using only sequence data and machine learning.

![#1589F0](https://placehold.it/10/1589F0/000000?text=+) **Note:** Notice that the model will classify an isolate as being *resistant* to rifampicin if at least one of the k-mers is absent in its genome. In fact, the rules capture the absence of the wild-type sequence, since all the variants at this locus were associated with resistance. Hence, to maximize the conciseness of the model, a single absence rule was used instead of using a presence rule for each variant.

For a detailed tutorial on model interpretation, see [here](./doc_interp.html).

## Predicting with the obtained model

You could use any k-mer counting tool, such as [Jellyfish](https://github.com/gmarcais/Jellyfish), [DSK](https://github.com/GATB/dsk), or [KMC](https://github.com/refresh-bio/KMC) to extract the k-mers present in a contig or read file and then, easily apply the model to new isolates.
This feature is currently not included in Kover, but will be added in future releases.
