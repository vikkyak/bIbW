ClusterNorm <-
  function(X,
           GeneFilter = NULL,
           Scaling = NULL,
           Method = NULL,
           Thresh = NULL) {
    if (is_sparse(X)) {
      X <- as(X, "dgCMatrix")
    }
    else {
      X <- as.matrix(X)
    }
    if (is.null(Scaling)) {
      Scaling <- colSums(X)
    }
    if (any(Scaling == 0)) {
      stop("cells should have non-zero library sizes or 'Scaling' values")
    }
    eigenV <- svd(X[GeneFilter, , drop = FALSE])[[1]]
    CountData <- X[GeneFilter, , drop = FALSE] / eigenV[1]
    ## Multifactor analysis of each conditions
    SF <- MultiFactor(CountData[, , drop = FALSE], Method, Thresh)
    Scaling <- Scaling / mean(eigenV)
    SF <- SF * sqrt(Scaling)
    X <- t(t(X[, , drop = FALSE]) / SF)
    list(X = X, SF = SF)
  }
