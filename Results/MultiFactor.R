MultiFactor <-
  function (Data, Method, Thresh, BPPARAM = SerialParam())
  {
    # mid quantile
    md <-
      apply(Data, 2, function(x)
        Qtools::midquantile(log(x[x > 0]), probs = c(0.25, 0.5, 0.75)))
    # A confidence interval for sample mid-quantiles can be obtained using confint.midquantile.
    
    #  lower bounds of the confidence intervals
    midL <-
      sapply(seq_along(md), function(x) {
        confint(md[[x]], level = 0.95)
      }[, "lower"])
    #  upper bounds of the confidence intervals
    midU <-
      sapply(seq_along(md), function(x) {
        confint(md[[x]], level = 0.95)
      }[, "upper"])
    
    rownames(midU)  <- rownames(midL)  <- c("25%", "50%" , "75%")
    
    mid <- list(midL, midU)
    
    mid <- lapply(seq_along(mid), function(x) {
      na.omit(mid[[x]])
    })
    
    mid <- lapply(seq_along(mid), function(x) {
      mid[[x]][apply(mid[[x]], 1, sd) != 0, , drop = FALSE]
    })
    
    plsc <- pls(t(mid[[1]]), t(mid[[2]]), ncomp = 1)
    pls1 <- plsc$variates$X
    pls2 <- plsc$variates$Y
    pc <- list(pls1,pls2)
    
    NonZeroGenes <- apply(Data, 1, function(x) {
      cc <- which(x > 0)
      list(rbind(cc, log(x[cc]), deparse.level = 0))
    })
    NonZeroGenes <- lapply(NonZeroGenes, function(x)
      matrix(unlist(x),
             nrow = 2))
    AllAvg <- unlist(lapply(NonZeroGenes, function(x)
      mean(x[2, ])))
    NonZeroCells <- apply(Data, 2, function(x)
      list(which(x >
                   0)))
    NonZeroCells <- lapply(NonZeroCells, function(x)
      c(unlist(x)))
    
    rm(Data, midL, midU, mid)
    gc()
    
    NonZeroGenesG <- list(NonZeroGenes)[rep(1, length(pc))]
    Bw <-
      bpmapply(
        FUN = BW,
        pc = pc,
        NonZeroGenes = NonZeroGenesG,
        Method,
        BPPARAM = BPPARAM
      )
    
    # Fuzzy weighted expression level of genes in the cells 
    # after upper and lower quantile evaluations
    Fw <- lapply(seq_along(pc), function(x) {
      LocAvg(NonZeroCells, NonZeroGenes, pc[[x]], Bw[, x],
             Thresh)
    })
    
    ## Individual
    LogScalingF1 <- rep(0, length(Fw[[1]]))
    for (j in 1:length(Fw[[1]])) {
      LogScalingF1[j] <- mean(Fw[[1]][[j]] - AllAvg[NonZeroCells[[j]]],
                             na.rm = T) - 1e-10
    }
    
    LogScalingF2 <- rep(0, length(Fw[[2]]))
    for (j in 1:length(Fw[[2]])) {
      LogScalingF2[j] <- mean(Fw[[2]][[j]] - AllAvg[NonZeroCells[[j]]],
                              na.rm = T) - 1e-10
    }
    
    (exp(LogScalingF1 - median(LogScalingF1)) + exp(LogScalingF2 - median(LogScalingF2)))
    # mean expression level of genes in the cells 
    # after upper and lower quantile evaluations
    # AdF <- list()
    # for (l in 1:length(Fw[[1]])) {
    #   AdF[[l]] <- (Fw[[1]][[l]] +  Fw[[2]][[l]]) / 2
    # }
    # # Log scaling factor
    # LogScalingF <- rep(0, length(AdF))
    # for (j in 1:length(AdF)) {
    #   LogScalingF[j] <- median(AdF[[j]] - AllAvg[NonZeroCells[[j]]],
    #                          na.rm = T) - 1e-10
    # }
    # # exp(LogScalingF - median(LogScalingF))
    # exp(LogScalingF - median(LogScalingF))
  }


