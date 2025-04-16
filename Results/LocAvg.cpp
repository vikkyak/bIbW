
#include <Rcpp.h>
using namespace Rcpp;

double WeightedAverage(NumericVector x, const double y, const double h, NumericVector z, const double a) {
  NumericVector d1 = (x - y)/h;
  NumericVector u = ifelse( abs(x-y) < a * h, exp(- 0.5 * d1 * d1), 0); /* u is the membership degree between the cells*/
  double sumu = sum(u);    /* sum of all u */
  double w = sum( u * z )/sumu;  /* weighted expression level of genes in the cells */
  return(w);
}

// [[Rcpp::export]]
List LocAvg(List NonZeroCell, List NonZeroGene, NumericVector pc, NumericVector b, const double Th) {
  int n = NonZeroCell.size();
  Rcpp::List output(n);
  
  for(int c=0; c<n; c++) {
    SEXP lc = NonZeroCell[c];
    NumericVector indexc(lc); indexc = indexc - 1;
    NumericVector outc(indexc.size());
    for(int g=0; g<indexc.size(); g++){
      SEXP lg = NonZeroGene[indexc[g]];
      NumericMatrix gM(lg);
      NumericVector indexg = gM(0,_) - 1;
      double bw = b[indexc[g]];
      outc[g] = WeightedAverage(pc[indexg], pc[c], bw, gM(1,_), Th);
    }
    output[c] = outc;
  }
  return(output);
}
