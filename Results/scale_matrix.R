scale_matrix <- function(matrix1, matrix_all) {
  indx1   <- which(as.matrix(matrix1) != 0)
  indxall <- which(as.matrix(matrix_all) != 0)
  matrix1[indx1]      <- log(matrix1[indx1])
  matrix_all[indxall] <- log(matrix_all[indxall])
  
  sc_fc <- rowSums(matrix1)/rowSums(matrix1 != 0) - rowSums(matrix_all)/rowSums(matrix_all != 0)
  exp(median(sc_fc, na.rm = TRUE))
  
}