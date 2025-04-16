
Nnaive <- function(counts, conditions = NULL) {
  #counts <- as.matrix(data); conditions <- c(rep(1,48), rep(2,44))
  G   <- nrow(counts)
  n   <- ncol(counts)
  cgs <- unique(conditions)
  
  
  SN <- sparseMatrix(G, n, x = 0)
  mg <- rep(NA, n)
  for(k in 1:length(cgs)) {
    ck <- conditions == cgs[k]
    cmatrix_k <- counts[, ck, drop = FALSE]
    sj_k <- colSums(cmatrix_k)
    mg[ck]   <- sj_k/median(sj_k)
    SN[, ck] <- t(t(cmatrix_k)/mg[ck])
  }
  
  if (length(cgs) <= 1) {
    reSN  <- SN
    remg2 <- mg
  } else {
    remg <- rep(NA, n)
    for(k in 1:length(cgs)) {
      ck <- conditions == cgs[k]
      remg[ck] <- mg[ck] * scale_matrix(SN[, ck], SN)
    }
    remg2 <- remg/median(remg)
    reSN  <- t(t(SN) * mg / remg2)
  }
  
  colnames(reSN) <- colnames(counts)
  rownames(reSN) <- rownames(counts)
  list(NormalizedData = reSN, scalingFactor = remg2)
  
}
