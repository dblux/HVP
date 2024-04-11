# RVP 

RVP (recursive variance partitioning) is a quantitative batch effect metric
that estimates the proportion of variance associated with batch effects in a
data set (the "RVP" value of a data set).

To determine whether batch effects are statistically significant in a data
set, a permutation test can be performed by setting `nperm` to the desired
number of permutations. We recommend performing at least 1000 permutations.

## Installation

### Install from Github

``` r
# Install package: devtools if not present
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("dblux/RVP")
```

## Usage

``` r
library(RVP)

# Specify number of samples to simulate for each batch-class group, with
# rows representing classes and columns representing batches.
crosstab <- matrix(5, 3, 2)

# Simulate microarray data with 100 features
simdata <- simulate_microarray(100, crosstab)

X <- simdata$X    # matrix with dimensions (nfeature, nsamples)
batch <- simdata$metadata$batch    # vector representing batch
class <- simdata$metadata$class    # vector representing class

res <- RVP(X, batch, class)
print(res$RVP)

# To perform permutation test
res <- RVP(X, batch, class, nperm = 1000)
print(res$p.value)
```

`RVP` is an S3 generic function; methods can be added for new
classes. S3 methods for class: `matrix`, `data.frame`, `SummarizedExperiment`,
`SingleCellExperiment` and `Seurat` are provided.

``` r
# "batch" and "class" are column names of colData in SingleCellExperiment object
res <- RVP(sce, "batch", "class", nperm = 1000)

# "batch" and "class" are column names of metadata in Seurat object
res <- RVP(seu, "batch", "class", nperm = 1000)
```
