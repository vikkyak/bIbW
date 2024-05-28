ClusterNorm <- function(X, GeneFilter = NULL, scaling = NULL, Method = NULL, cutoff = NULL ){
  
  if (is_sparse(X)) {
    X <- as(X, "dgCMatrix")
  }
  else {
    X <- as.matrix(X)
  }
  if (is.null(scaling)) {
    scaling <- colSums(X)
  }
  if (any(scaling == 0)) {
    stop("cells should have non-zero library sizes or 'scaling' values")
  }
  eigenV <- svd(X[GeneFilter, , drop = FALSE])[[1]]
  CountData <- X[GeneFilter, , drop = FALSE] / eigenV[1]
  SF <- MultiFactor(CountData[, , drop = FALSE], Method, cutoff)
  scaling <- scaling / mean(eigenV)
  SF <- SF * sqrt(scaling)
  X <- t(t(X[, , drop = FALSE]) / SF)
  list(X = X, SF = SF)
}
