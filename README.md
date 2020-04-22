#normR - normR obeys regime mixture rules
## Normalization and Difference Calling for Next Generation Sequencing (NGS)
Experiments via Joint Multinomial Modeling
---

Two NGS tracks are modeled simultaneously by fitting a binomial mixture model
on mapped read counts.  In the first counting process, a desired smoothing
kernel (bin size) and read characteristic threshold (quality, SAMFLAG) can be
specified.  In a second step a binomial mixture model with a user-specified
number of components is fit to the data.  The fit yields different enrichment
regimes in the supplied NGS tracks.  Log-space computation is done in C/C++
where [OpenMP](http://openmp.org) enables for fast parallel computation.


### Release Version

The master branch is always in sync with the [normR Bioconductor release](
http://bioconductor.org/packages/devel/bioc/html/normr.html) and the [normR
github Bioconductor mirror]( https://github.com/Bioconductor-mirror/normr).  A
R 3.2 compliant version can be found in the [normR R3.2 tree](
https://github.com/your-highness/normR/tree/R3.2).


### Installation

To install normR from the release repository, easiest way is to use
Bioconductor or devtools:

```R
#install dependencies
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("bamsignals", suppressUpdates=T)
#fetch current normR version from github
install.packages("devtools")
require(devtools)
devtools::install_github("your-highness/normr")
```

### Usage

See the
[vignette](
https://cdn.rawgit.com/your-highness/normR/development/inst/doc/normr.html) 
for a toy example on normR usage. The documentation of routines can be accessed
  from with R with ``?``.

#### Use cases

* ChIP-seq normalization / enrichment calling with an Input experiment (Whole
  Cell Extract, H3/IgG ChIP-seq)

* ChIP-seq differential enrichment calling for two different antigens in same
  sample population

* ChIP-seq identification of enrichment regimes to investigate on sample
  heterogeneity

* RNA-seq differential expression calling

* ChIP-seq differential enrichment calling in two different samples (be aware
  of CNVs!)

* CNV identification


#### Useful links

Be sure to check out the following amazing github projects for your upcoming
NGS magic:

[bamsignals](https://github.com/lamortenera/bamsignals) - Efficient Counting in
Indexed Bam Files for Single End and Paired End NGS Data

[EpicSeg](https://github.com/lamortenera/epicseg) - Chromatin Segmentation
Based on a Probabilistic Multinomial Model for Read Counts

[kfoots](https://github.com/lamortenera/kfoots) - Fit Multivariate Discrete
Probability Distributions to Count Data

[deepTools](https://github.com/fidelram/deepTools) - User-Friendly Tools for
Normalization and Visualization of Deep-Sequencing Data
