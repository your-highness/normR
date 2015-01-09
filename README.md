diffr
=====

## Normalization and difference calling for Next Generation Sequencing (NGS) experiments via joint multinomial modeling

Functions for normalization and difference calling in NGS experiment setups. A
binomial mixture model with a given number of components is fit and used for
identifying enriched or depleted regions in two given data tracks. Log-space
multinomial model is fit by Expectation maximization in C/C++.

The background component is assumed to be the component with lowest mean and is
used to compute treatment fold change and P-Values for statistical significance
of enrichment.

## Installation

To install diffr from the working repository, easiest is using devtools:
```R
install.packages("devtools")
devtools::install_github("your-highness/diffr")
```
## Use cases

RNA-seq differential expression calling

ChIP-seq normalization with Input experiment, ChIP-seq differential enrichment
calling for different proteins in same cell type, ChIP-seq differential
enrichment calling for same proteins in different cell types, ChIP-seq
differential enrichment calling for different proteins in different cell types

ChIP-seq identification of enrichment regimes to investigate on two vs one
histone tail modifications or cell population mixture effects
