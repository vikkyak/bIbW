BW <- function( NonZeroGenes, pc, Method = NULL){ 
  
  
  if (Method == "KernSmooth") {
  hh <- sapply( NonZeroGenes, function(x) {
    KernSmooth::dpik(pc[x[1,]])
  })
}
  else if (Method == "SJ")  {
    hh <- sapply( NonZeroGenes, function(x) {
      h <- tryCatch({
        bw.SJ(pc[x[1, ]], method = "ste")
      }, error = function(e) {
        return(bw.nrd0(pc[x[1, ]]))
      })
    })
  }
  
  else if (Method == "ucv") {
    hh <- sapply( NonZeroGenes, function(x) {
      h <- tryCatch({
        bw.ucv(pc[x[1, ]])
      }, error = function(e) {
        return(bw.nrd0(pc[x[1, ]]))
      })
    })
  }
  
  else if (Method == "bcv") {
    hh <- sapply( NonZeroGenes, function(x) {
      h <- tryCatch({
        bw.bcv(pc[x[1, ]])
      }, error = function(e) {
        return(bw.nrd0(pc[x[1, ]]))
      })
    })
  }
  
  else {
    
    hh <- sapply( NonZeroGenes, function(x) {
      h <- tryCatch({
        bw.nrd(pc[x[1, ]])
      }, error = function(e) {
        return(bw.nrd0(pc[x[1, ]]))
      })
    })
  }
  hh
}