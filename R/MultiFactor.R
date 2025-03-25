MultiFactor <-function (Data, Method, cutoff, BPPARAM = SerialParam())
{
  
  # mid quantile
  md <-
    apply(Data, 2, function(x)
      Qtools::midquantile(log(x[x > 0]), probs = c(0.25, 0.5, 0.75)))
  
  midL <-
    sapply(seq_along(md), function(x) {
      confint(md[[x]], level = 0.95)
    }[, "lower"])
  
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
  
  pca <- lapply(seq_along(mid), function(x) {
    stats::prcomp(t(mid[[x]]), scale = TRUE)
  })
  
 pc <- lapply(seq_along(pca), function(x) {
    pca[[x]]$x[, 1]
  })
  
  NonZeroGenes <- apply(Data, 1, function(x) {
    cc <- which(x > 0)
    list(rbind(cc, log(x[cc]), deparse.level = 0))
  })
  NonZeroGenes <- lapply(NonZeroGenes, function(x)
    matrix(unlist(x),
           nrow = 2))
  AllAvg <- unlist(lapply(NonZeroGenes, function(x)
    mean(x[2,])))
  NonZeroCells <- apply(Data, 2, function(x)
    list(which(x >
                 0)))
  NonZeroCells <- lapply(NonZeroCells, function(x)
    c(unlist(x)))
  
  rm(Data, midL, midU, mid, pca)
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
  
  ## Fuzzy weights
  Fw <- lapply(seq_along(pc), function(x) {
    LocAvg(NonZeroCells, NonZeroGenes,pc[[x]], Bw[, x],
                    cutoff)
  })
  
  ## Average Defuzzification
  AdF <- list()
  for (l in 1:length(Fw[[1]])) {
    AdF[[l]] <- (Fw[[1]][[l]] +  Fw[[2]][[l]]) / 2
  }
  
  LogScalingF <- rep(0, length(AdF))
  for (j in 1:length(AdF)) {
    LogScalingF[j] <- mean(AdF[[j]] - AllAvg[NonZeroCells[[j]]],
                           na.rm = T) - 1e-10
  }
  exp(LogScalingF - median(LogScalingF))
  
  }
