#' @useDynLib abcGWAS
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics plot
#' @importFrom stats dist ecdf quantile rgamma
#' @importFrom utils head read.table setTxtProgressBar txtProgressBar write.table setTxtProgressBar txtProgressBar 
#' @importFrom ggplot2 aes aes_string facet_wrap geom_histogram geom_path geom_vline ggplot ggtitle theme 
NULL
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
