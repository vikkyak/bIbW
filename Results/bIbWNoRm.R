bIbWNoRm <-
  function (Data,
            Clusters = NULL,
            Scaling = NULL,
            FilterExpression = 0,
            FilterCellNum = 10,
            Method = c("KernSmooth", "SJ", "ucv", "bcv", "nrd"),
            Thresh = 1.5,
            MF = NULL,
            BPPARAM = SerialParam())
  {
    Method <- match.arg(Method)
    nCells <- ncol(Data)
    col.order <- colnames(Data)
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
    if (!is.null(Scaling) && length(Scaling) != ncol(Data)) {
      stop("'length(Scaling)' should be equal to 'ncol(x)'")
    }
    
    GroupData <- GroupScale <- vector("list", length(indices))
    
    Levels <- unique(Clusters)
    
    GroupData <- lapply(seq_along(Levels), function(x) {
      Data[, which(Clusters == Levels[x])]
    })
    if (length(GroupData)==1) {
      plsData <- GroupData
    }
    else{
    a <- combn(seq_along(GroupData), 2)
    corrData <- lapply(1:ncol(a), function(x) {
      ind1 <- a[1, x]
      ind2 <- a[2, x]
      respls <- pls(GroupData[[ind1]], GroupData[[ind2]], ncomp = 1)
      corr <- cor(respls$variates$X, respls$variates$Y)
      # plsdata <- as.vector(corr) * (cbind(GroupData[[ind1]], GroupData[[ind2]]))
    })
    
    # maxcor <- round(max(unlist(corrData)),3)
    maxcor <- round(mean(unlist(corrData)),3)
    plsData <- lapply(seq_along(Levels), function(x) {
      D <-maxcor*GroupData[[x]]
    })
    }
    
    
    
    GroupScale <- lapply(seq_along(Levels), function(x) {
     Scaling[x]
    })
    
    Genes <- rownames(Data)
    SeqDepthGrouped <- lapply(seq_along(Levels), function(x) {
      colSums(Data[, which(Clusters == Levels[x])])
    })
    
    NumZerosCellGrouped <- lapply(seq_along(Levels), function(x) {
      colSums(plsData[[x]] != 0)
    })
    if ((sum(do.call(c, SeqDepthGrouped) == 10000) / (nrow(Data) *
                                                      ncol(Data))) >= 0.8) {
      warning(
        "More than 80% of data is zeros.  \n        Check data (remove low quality cells prior). \n
            You may also need to adjust the filtering criteria \n        parameters FilterExpression and FilterCellNum."
      )
    }
    if (any(do.call(c, NumZerosCellGrouped) <= 100)) {
      warning(
        "At least one cell/sample has less than 100 genes detected (non-zero). \n
            Check data or filtering criteria. \n
            bIbW may not be appropriate for your data."
      )
    }
    message("Gene filtering is performed within each condition or cluster")
    GeneZerosGropued <- lapply(seq_along(Levels), function(x) {
      rowSums(plsData[[x]] != 0)
    })
    MedExprGrouped <- lapply(seq_along(Levels), function(x) {
      apply(plsData[[x]], 1, function(c)
        median(c[c != 0]))
    })
    GeneFilterGrouped <- lapply(seq_along(Levels), function(x) {
      names(which(
        GeneZerosGropued[[x]] >= FilterCellNum & MedExprGrouped[[x]] >=
          FilterExpression
      ))
    })
    checkGeneFilter <- vapply(seq_along(Levels), function(x) {
      length(GeneFilterGrouped[[x]])
    }, FUN.VALUE = numeric(1))
    if (any(checkGeneFilter < 100)) {
      stop(
        "At least one Cluster has less then 100 genes that pass the specified filter.
  Check  your data or filtering criteria. \n bIbW may not be appropriate for your data."
      )
    }
    GeneFilterOUT <- lapply(seq_along(Levels), function(x) {
      names(which(
        GeneZerosGropued[[x]] < FilterCellNum | MedExprGrouped[[x]] <
          FilterExpression
      ))
    })
    names(GeneFilterOUT) <- paste0("GenesFilteredOutGroup",
                                   unique(Clusters))
    NM <- lapply(seq_along(Levels), function(x) {
      message(
        paste0(
          length(GeneFilterOUT[[x]]),
          " genes in Cluster ",
          Levels[x],
          " will not be included in the normalization due to \n the defined filtering criteria."
        )
      )
    })
    
    
    
    if (is.null(Clusters) | length(unique(Clusters)) ==
        1) {
      SF <- MultiFactor(Data[GeneFilterGrouped[[1]],], Method, Thresh)
      Data <- as(t(t(Data) / SF), "sparseMatrix")
      ScaleFactors <- SF
    }
    else {
      BiGroupedData <-
        bpmapply(
          FUN = ClusterNorm,
          X = plsData,
          GeneFilter = GeneFilterGrouped,
          Scaling = GroupScale,
          SIMPLIFY = FALSE,
          USE.NAMES = FALSE,
          Method = Method,
          Thresh = Thresh,
          MF = MF,
          BPPARAM = BPPARAM
        )
      
      
      Data <- lapply(BiGroupedData, "[[", i = "X")
      ScaleFactors <- lapply(BiGroupedData, "[[", i = "SF")
      len <- length(BiGroupedData)
      
      SF  <- unlist(lapply(seq_len(len), function(x) {
        BiGroupedData[[x]]$SF}))[col.order]
      
      NormData <- do.call(cbind, lapply(seq_len(len), function(x) {
        cbind(BiGroupedData[[x]]$X)
      }))
      NormData <- NormData[, col.order]
      
      Index <- lapply(seq_len(len), function(x) {
        BiGroupedData[[x]]$X != 0
      })
      
      SparseCount <-
        lapply(seq_len(len), function(x) {
          as(BiGroupedData[[x]]$X, "sparseMatrix")
        })
      
      for (l in 1:len) {
        SparseCount[[l]]@x = log(SparseCount[[l]]@x)
      }
      
     Ref <- lapply(seq_len(len), function(x) {
        rowSums(SparseCount[[x]][, , drop = FALSE]) / rowSums(Index[[x]])
      })
      
     Ref <-  Reduce("+",Ref)
      
      ScaleFactors <- lapply(seq_len(len), function(x) {
        (BiGroupedData[[x]]$SF) * exp(median(
          rowSums(SparseCount[[x]][, , drop = FALSE]) /
            rowSums(Index[[x]][, , drop = FALSE]) -Ref,
          na.rm = T
        ))
      })
      
      # ScaleFactors <- sapply(seq_len(len), function(x) {
      #   (BiGroupedData[[x]]$SF) * exp(median(
      #     rowSums(SparseCount[[x]][, , drop = FALSE]) /
      #       rowSums(Index[[x]][, , drop = FALSE]) -Ref,
      #     na.rm = T
      #   ))
      # })
      
      
      ScaleFactors <- unlist(ScaleFactors)[col.order]
      ScaleFactors <- ScaleFactors / median(ScaleFactors)
      # ScaleFactors <- ScaleFactors / geometric.mean(ScaleFactors)
      
      Data <- as(t(t(NormData) * SF/ScaleFactors), "sparseMatrix")
  
      # Data <-
      #   qsmooth::qsmooth(object = as.matrix(NormData[, , drop = FALSE]) , group_factor = Clusters)
    }
    list(NormalizedData = Data, ScaleFactors = ScaleFactors)
  }
