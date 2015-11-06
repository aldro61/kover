---
title: Getting started with Kover
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

This genomes and the metadata were obtained from: Merker, Matthias, et al. "Evolutionary history and global spread of the Mycobacterium tuberculosis Beijing lineage." *Nature genetics* 47.3 (2015): 242-249.

## Creating a dataset

## Learning a model