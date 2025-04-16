filter.genes <- function(RawCounts){
  raw.counts<-RawCounts
  #remove genes that are never expressed
  keep<-(rowSums(raw.counts)>3)
  raw.counts<-raw.counts[keep,]
  Rpm<-apply(raw.counts, 2, function(x) 1e6*(x/sum(x))) #rpm(RawCounts)$counts
  Mean<-rowMeans(Rpm)
  Select<-names(Mean)[which(Mean>median(Mean))]
  return(list(Select=Select, Mean=Mean))
}