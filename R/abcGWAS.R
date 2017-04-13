#' @useDynLib abcGWAS
#' @importFrom Rcpp sourceCpp
NULL
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
