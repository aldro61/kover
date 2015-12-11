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

Kover is an *out-of-core* implementation of the Set Covering Machine learning algorithm.
Given two sets of whole genomes, it attempts to find a simple model, that relies on a small number of short DNA sequences ([k-mers](https://en.wikipedia.org/wiki/K-mer)), that accurately discriminate them.
Kover can be used to obtain interpretable models of phenotypes.


## Survey of features

We provide a command line interface for Kover. Kover's features are regrouped in two tools:

* [Kover dataset tools](doc_dataset.html): Manipulate and prepare genomic data for Kover
* [Kover learning engine](doc_learning.html): Apply the machine learning algorithm to learn models

The back-end of Kover is bundled as a Python package. All the actions available in the command line interface can be
performed using this package. This is not officially supported, but basic documentation is provided in the 
docstring of some functions.


## Licence

Kover is open-source software released under the [GPLv3 licence](http://www.gnu.org/licenses/gpl-3.0.html).


## Getting started

To get started, see these two topics:

1. {{site.data.urls.doc_installation.link}}
2. {{site.data.urls.doc_example.link}}