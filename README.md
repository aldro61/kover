<img src="./kover.png" height="50" />2.0

[![DOI](https://zenodo.org/badge/20289/aldro61/kover.svg)](https://zenodo.org/badge/latestdoi/20289/aldro61/kover)
[![Build Status](https://travis-ci.org/aldro61/kover.svg?branch=kover2)](https://travis-ci.org/aldro61/kover)

Kover is an *out-of-core* implementation of rule-based machine learning algorithms that has been tailored for genomic biomarker discovery. It produces highly interpretable models, based on k-mers, that explicitly highlight genotype-to-phenotype associations.

## Introduction

Understanding the relationship between the genome of a cell and its phenotype is a central problem in precision medicine. 
Nonetheless, genotype-to-phenotype prediction comes with great challenges for machine learning algorithms that limit their use in this setting. The high dimensionality of the data tends to hinder generalization and challenges the scalability of most learning algorithms. Additionally, most algorithms produce models that are complex and difficult to interpret. We alleviate these limitations by proposing strong performance guarantees, based on sample compression theory, for rule-based learning algorithms that produce highly interpretable models. We show that these guarantees can be leveraged to accelerate learning and improve model interpretability. Our approach is validated through an application to the genomic prediction of antimicrobial resistance, an important public health concern. Highly accurate models were obtained for 12 species and 56 antibiotics, and their interpretation revealed known resistance mechanisms, as well as some potential new ones. An open-source disk-based implementation that is both memory and computationally efficient is included with this work. The implementation is turnkey, requires no prior knowledge of machine learning, and is complemented by comprehensive tutorials.

> Drouin, A., Letarte, G., Raymond, F., Marchand, M., Corbeil, J., & Laviolette, F. (2019). Interpretable genotype-to-phenotype classifiers with performance guarantees. Scientific Reports, 9(1), 4071. [[PDF]](https://www.nature.com/articles/s41598-019-40561-2)

> Drouin, A., Giguère, S., Déraspe, M., Marchand, M., Tyers, M., Loo, V. G., Bourgault, A. M., Laviolette, F. & Corbeil, J. (2016). Predictive computational phenotyping and biomarker discovery using reference-free genome comparisons. BMC Genomics, 17(1), 754. [[PDF]](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2889-6)


### Video lecture:
The Set Covering Machine implementation in Kover was featured in the following video lecture:

Interpretable Models of Antibiotic Resistance with the Set Covering Machine Algorithm, Google, Cambridge, Massachusetts (February 2017)

[![Google tech talk](https://img.youtube.com/vi/uomMdBdEwnk/0.jpg)](https://www.youtube.com/watch?v=uomMdBdEwnk)

## Installation

You can use either of the following options:
* Docker image with Kover preinstalled (https://hub.docker.com/r/aldro61/kover)
* Manual installation: http://aldro61.github.io/kover/doc_installation.html

## Tutorials

For tutorials on how to use Kover with your data, see: http://aldro61.github.io/kover/doc_tutorials.html

## Documentation

The documentation can be found at: http://aldro61.github.io/kover/

## Contact

If you need help using Kover, please use [Biostars](https://www.biostars.org/p/194520/). To report a bug, please create an issue on GitHub.
