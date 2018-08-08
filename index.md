---
title: Introduction
tags: 
  - "getting-started"
type: first_page
toc: false
homepage: true
published: true
---

## Overview 

Kover is an out-of-core implementation of rule-based machine learning algorithms that has been tailored for genomic biomarker discovery. It produces highly interpretable models, based on k-mers, that explicitly highlight genotype-to-phenotype associations.

Given two groups of phenotipically distinct individuals represented by their genomes, Kover seeks an intelligible model that
accurately discriminates them. The obtained models are arrangements of rules that capture
the presence or absence of [k-mers](https://en.wikipedia.org/wiki/K-mer).


For example, when applied to 462 [*C. difficile*](https://en.wikipedia.org/wiki/Clostridium_difficile_(bacteria)) isolates divided into
two groups: resistant or sensitive to [Azithromycin](https://en.wikipedia.org/wiki/Azithromycin), Kover found that the
following model is a good predictor of resistance to this drug:

```
Presence(AGCCAGGTTCTTCATTTAAGATGCTAACTTC)
OR
Presence(CTTAAGCTGCCAGCGGAATGCTTTCATCCTA)
OR
Presence(AAGTCGCCCTTTTTTAAGGATACGGCGGTAT)
```

## Survey of features

Kover currently implements two learning algorithms: 1) Set Covering Machines, 2) Classification Trees. Both algorithms can be trained out-of-core, i.e., without loading the entire dataset into the computer's memory. The implementation are fast and rely on data compression and atomic CPU instructions to accelerate computation.

We provide tight performance guarantees on the obtained models and use them to accelerate learning through a model selection strategy called [bound selection](./doc_learning.html/#risk-bound-selection).

A command line interface is provided. It consists in two main modules, [kover dataset](doc_dataset.html) and [kover learn](doc_learning.html). Kover dataset provides data manipulation utilities and kover learn is an interface to the machine learning algorithm.

## Licence

Kover is open-source software released under the [GPLv3 licence](http://www.gnu.org/licenses/gpl-3.0.html).


## Getting started

To get started, see these two topics:

1. {{site.data.urls.doc_installation.link}}
2. {{site.data.urls.doc_example.link}}