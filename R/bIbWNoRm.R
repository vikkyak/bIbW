
#' Within and Between (bIbW) Samples Normalization of RNA-seq data
#' @export
#' @param Data as a Matrix with rows as genes (features) and column as cells (samples).
#' @param Conditions and Methods for variance selection (Bandwidth).
#' @returns Output Normalized data and scale factor.
#' @examples
#' library(scran)
#' library(qsmooth)
#' library(Qtools)
#' library(dynutils)
#' library(BiocParallel)
#' set.seed(12345)
#' G <- 2000; n <- 600 # G: number of genes, n: number of cells
#' mu <- rgamma(G, shape = 2, rate = 2)
#' NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu)
#' SimulatedData <- sapply(1:n, NB_cell)
#' colnames(SimulatedData) <- paste("c", 1:n, sep = "_")
#' rownames(SimulatedData) <- paste("g", 1:G, sep = "_")
#' Conditions = rep(c(1,2), each= 300)
#' Clusters <- quickCluster(SimulatedData)
#' Result <- bIbWNoRm(Data = SimulatedData, Clusters = Clusters, Method= "KernSmooth")


bIbWNoRm <- function (Data,  Clusters = NULL, MaxClusterSize = 3000, scaling = NULL, filter = FALSE, FilterExpression = 0, 
                      FilterCellNum = 10, numGeneforEst = 2000, divideforFast = TRUE, 
                      numDivide = NULL, Method = c("KernSmooth","SJ","ucv","bcv", "nrd" ), cutoff = 1.5, BPPARAM = SerialParam()) 
{
 
 
   Method <- match.arg(Method)
  nCells <- ncol(Data)
  if (is.null(Clusters)) {
    Clusters <- integer(nCells)
  }
  
  if (nCells != length(Clusters)) {
    stop("'ncol(Data)' is not equal to 'length(Clusters)'")
  }
  
  indices <- split(seq_along(Clusters), Clusters)
  
  if (length(indices) == 0L || any(lengths(indices) == 0L)) {
    stop("zero cells in one of the clusters")
  }
  if (!is.null(scaling) && length(scaling) != ncol(Data)) {
    stop("'length(scaling)' should be equal to 'ncol(x)'")
  }
  
  GroupData <- GroupScale <- vector("list", length(indices))
  
  Levels <- unique(Clusters)
  
  GroupData <- lapply(seq_along(Levels), function(x) {
    Data[, which(Clusters == Levels[x])]
  })
  GroupScale <- lapply(seq_along(Levels), function(x) {
    scaling[x]
  })
  
  Genes <- rownames(Data)
  SeqDepthGrouped <- lapply(seq_along(Levels), function(x) {
    colSums(Data[, which(Clusters == Levels[x])])
  })
  
  NumZerosCellGrouped <- lapply(seq_along(Levels), function(x) {
    colSums(GroupData[[x]] != 0)
  })
  if ((sum(do.call(c, SeqDepthGrouped) == 10000)/(nrow(Data) * 
                                                  ncol(Data))) >= 0.8) {
    warning("More than 80% of your data is zeros.  \n        Check the quality of your data (remove low quality cells prior to running SCnorm). \n        You may need to adjust the filtering criteria for SCnorm using\n        parameters FilterExpression and FilterCellNum. \n        It could also be the case that SCnorm is not be appropriate for your data (see vignette for details).")
  }
  if (any(do.call(c, NumZerosCellGrouped) <= 100)) {
    warning("At least one cell/sample has less than 100 genes detected (non-zero). \n        Check the quality of your data or filtering criteria. \n        SCnorm may not be appropriate for your data (see vignette for details).")
  }
  message("Gene filter is applied within each condition.")
  GeneZerosGropued <- lapply(seq_along(Levels), function(x) {
    rowSums(GroupData[[x]] != 0)
  })
  MedExprGrouped <- lapply(seq_along(Levels), function(x) {
    apply(GroupData[[x]], 1, function(c) median(c[c != 0]))
  })
  GeneFilterGrouped <- lapply(seq_along(Levels), function(x) {
    names(which(GeneZerosGropued[[x]] >= FilterCellNum & MedExprGrouped[[x]] >= 
                  FilterExpression))
  })
  checkGeneFilter <- vapply(seq_along(Levels), function(x) {
    length(GeneFilterGrouped[[x]])
  }, FUN.VALUE = numeric(1))
  if (any(checkGeneFilter < 100)) {
    stop("At least one Cluster has less then 100 genes that pass the specified filter. Check the quality of your data or filtering criteria. \n       SCnorm may not be appropriate for your data (see vignette for details).")
  }
  GeneFilterOUT <- lapply(seq_along(Levels), function(x) {
    names(which(GeneZerosGropued[[x]] < FilterCellNum | MedExprGrouped[[x]] < 
                  FilterExpression))
  })
  names(GeneFilterOUT) <- paste0("GenesFilteredOutGroup", 
                                 unique(Clusters))
  NM <- lapply(seq_along(Levels), function(x) {
    message(paste0(length(GeneFilterOUT[[x]]), " genes in Cluster ", 
                   Levels[x], " will not be included in the normalization due to \n             the specified filter criteria."))
  })
  message("A list of these genes can be accessed in output, \n    see vignette for example.")
  
  
  if (is.null(Clusters) | length(unique(Clusters)) == 
      1) {
    SF <- MultiFactor(Data[GeneFilterGrouped[[1]], ], Method, cutoff)
    Data <- t(t(Data)/SF)
  }
  else {
    
    
    BiGroupedData <- bpmapply(FUN =ClusterNorm, X = GroupData, GeneFilter = GeneFilterGrouped,
                              scaling = GroupScale, SIMPLIFY=FALSE, USE.NAMES=FALSE, Method = Method, cutoff = 2, BPPARAM = BPPARAM)
    
    
    Data <- lapply( BiGroupedData, "[[", i = "X")
    ScaleFactors <- lapply(BiGroupedData, "[[", i = "SF")
    
    len <- length(BiGroupedData)  
    NormData <- do.call(cbind, lapply(seq_len(len), function(x) {
      cbind(BiGroupedData[[x]]$X)
    }))
    
    indxall <- lapply(seq_len(len), function(x) {
      BiGroupedData[[x]]$X !=0})
    
    SparseCount <- lapply(seq_len(len), function(x) {as(BiGroupedData[[x]]$X, "sparseMatrix")})
    
    for (l in 1: len) {
      SparseCount[[l]]@x = log(SparseCount[[l]]@x)
    }
    
    refall <- lapply(seq_len(len), function(x) {
      Matrix::rowSums(SparseCount[[x]][, , drop = FALSE])/rowSums(indxall[[x]])})
    
    refall <-  Reduce("+", refall)
    
    ScaleFactors <- sapply(seq_len(len), function(x) {
      (BiGroupedData[[x]]$SF) * exp(median(Matrix::rowSums(SparseCount[[x]][, , drop = FALSE]) /
                                             rowSums(indxall[[x]][, , drop = FALSE]) - refall, na.rm = T))
    })
    
    ScaleFactors <- unlist(ScaleFactors)
 
    ScaleFactors <-ScaleFactors/median(ScaleFactors)
    
    Data <- qsmooth::qsmooth(object = as.matrix(NormData[, , drop = FALSE]) , group_factor = Clusters)
  }
  list(NormalizedData = Data, ScaleFactors = ScaleFactors)
}
