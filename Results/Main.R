rm(list = ls())
library(BiocParallel)
library(scran)
library(Rcpp)
library(dynutils)
library(psych)
library(mixOmics)
source("~/Desktop/NormWithPLS/NewCor/SourceFun.R")
set.seed(12345)
G <- 2000; n <- 600 # G: number of genes, n: number of cells
mu <- rgamma(G, shape = 2, rate = 2)
NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu)
SimulatedData <- sapply(1:n, NB_cell)
colnames(SimulatedData) <- paste("c", 1:n, sep = "_")
rownames(SimulatedData) <- paste("g", 1:G, sep = "_")
Conditions = rep(c(1,2), each= 300)
Clusters <- quickCluster(SimulatedData)
bIbW <- bIbWNoRm(Data = SimulatedData, Clusters = Clusters, Method= "KernSmooth", MF = FALSE)




