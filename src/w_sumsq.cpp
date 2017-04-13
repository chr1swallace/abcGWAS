#include <Rcpp.h>
using namespace Rcpp;

/* Compute the sum of squares between x and y, weighted by w */

// [[Rcpp::export]]
double wsumsq( const NumericVector& xx,  const NumericVector& yy,  const NumericVector& ww){
  int k;
  double sumw =0.0;
  double sumxy =0.0;
  double sumsq =0.0;
  	/*  get dimensions of data */
  const int n = xx.size();
  /* for loop to calculate sums, the for loop is needed as we have a check for pairwise missing: ISNAN  */
  for (k = 0; k < n; k++){
    if (!ISNAN(xx[k]) && !ISNAN(yy[k])){
      sumw += ww[k];
      sumxy += ww[k]*(xx[k]-yy[k])*(xx[k]-yy[k]);
      // Rcpp::Rcout << xx[k] << ' ' << yy[k] << ' ' << sumxy << ' ' << sumw << std::endl;
    }
  }
  /*  final correlation */
  sumsq = sumxy/sumw;
  return sumsq;
}

// [[Rcpp::export]]
NumericVector wsumsqmat( const NumericVector& xx,  const NumericMatrix& yy,  const NumericVector& ww){
  /*  get dimensions of data */
  const int n = xx.size();
  const int m = yy.ncol();
  // Rcpp::Rcout << n << ' ' << m << std::endl;
  int i,j;
  NumericVector sumsq(m);
  /* for loop to calculate sums, the for loop is needed as we have a check for pairwise missing: ISNAN  */
  for (j = 0; j < m; j++){
    double sumw =0.0;
    double sumxy =0.0;
    for(i=0; i<n; i++) {
      if (!ISNAN(xx[i]) && !ISNAN(yy(i,j))){
	sumw += ww[i];
	sumxy += ww[i]*(xx[i]-yy(i,j))*(xx[i]-yy(i,j));
      }
    }
    /*  final correlation */
  // Rcpp::Rcout << sumxy << ' ' << sumw << std::endl;
    sumsq[j] = sumxy/sumw;
  }
  return sumsq;
}
