<img src="http://graal.ift.ulaval.ca/adrouin/kover.png" height="50" />

Kover is an *out-of-core* implementation of the Set Covering Machine learning algorithm. It compares large sets of wholes genomes and produces highly interpretable models of phenotypes.

# Introduction

Drouin, A., Giguère, S., Déraspe, M., Marchand, M., Tyers, M., Loo, V. G., Bourgault, A. M., Laviolette, F. & Corbeil, J. (2015). Predictive computational phenotyping and biomarker discovery using reference-free genome comparisons. submitted.
> Case-control studies compare groups of related genomes with the aim of identifying biomarkers that are predictive of a phenotype.
> Recent advances in next-generation sequencing have led to a tremendous increase in the scale of such studies.
> This trend will persist, motivating the need for efficient computational biomarker discovery and validation methodologies.
> We present a novel reference-free method for genomic biomarker discovery that produces uncharacteristically sparse models.
> It relies on a k-mer representation of genomes and on a greedy machine learning algorithm.
> The method successfully predicted the antibiotic resistance of four common human pathogens: *C. difficile*, *M. tuberculosis*, *P. aeruginosa* and *S. pneumoniae*, and yielded computational models for 17 commonly used antibiotics.
> We show that the models are accurate, faithful to the biological pathways targeted by the antibiotics and that they provide insight into the process of resistance acquisition. 
> Kover, an efficient implementation of our method, can readily scale to large genomic datasets.
> It is open-source and can be obtained from http://github.com/aldro61/kover.

Until the journal paper (more detailed) is published, please refer to this introductory paper:

Drouin, A., Giguère, S., Déraspe, M., Laviolette, F., Marchand, M., & Corbeil, J. (2015). Greedy Biomarker Discovery in the Genome with Applications to Antimicrobial Resistance. arXiv preprint arXiv:1505.06249. http://arxiv.org/pdf/1505.06249v1.pdf

# Installation

For installation instructions, see: http://aldro61.github.io/kover/doc_installation.html

# Example

For an example of use, see: http://aldro61.github.io/kover/doc_example.html

# Documentation

The documentation can be found at: http://aldro61.github.io/kover/
