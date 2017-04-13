#' @importFrom snpStatsWriter write.plink
#' @param XX SnpMatrix
#' @export
ldak.weights <- function(XX) {
    d <- tempdir()
    write.plink(file.base=file.path(d,"ldak"), snps=XX, phenotype=rep(1,nrow(XX)),chromosome=1)
    c1 <- paste("/scratch/wallace/abc/bin/run-ldak.sh",d)
    system(c1)
    weight_table <- read.table(file.path(d,"weightsALL"))
    weights <- weight_table$V1[unlist(lapply(snpswithp,function(x){which(weight_table$V5 == x)}))]
    names(weights) <- colnames(XX)
    weights
}

