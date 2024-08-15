# bIbW: Normalization of Single-cell RNA-seq Data using partial least squares (PLS) with Adaptive Fuzzy Weight

## Installation
``` r
library(devtools)
install_github("vikkyak/bIbW")
```
## Usage
``` r
library(bIbW)
```
### A simulated count matrix generated from negative binomial distributions

``` r
set.seed(123)
library(scran)
library(BiocParallel)
library(Qtools)
library(dynutils)
library(qsmooth)
set.seed(12345)
G <- 2000; n <- 600 # G: number of genes, n: number of cells
mu <- rgamma(G, shape = 2, rate = 2)
NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu)
SimulatedData <- sapply(1:n, NB_cell)
colnames(SimulatedData) <- paste("c", 1:n, sep = "_")
rownames(SimulatedData) <- paste("g", 1:G, sep = "_")
Conditions = rep(c(1,2), each= 300)
Clusters <- quickCluster(SimulatedData)
Result <- bIbWNoRm(Data = SimulatedData, Clusters = Clusters, Method= "KernSmooth")
```


