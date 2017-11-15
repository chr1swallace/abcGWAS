#' Calculate SNP weights from a SnpMatrix object
#'
#' This function calls the external executable ldak, which is assumed
#' to exist in your PATH environment variable.
#'
#' LDAK can be downloaded from http://dougspeed.com/ldak/.
#'
#' We used version 4.9.  Later versions may not work if the names of
#' core arguments changed, but previous versions appear to be kept on
#' the ldak website.
#' @importFrom snpStats write.plink
#' @param XX SnpMatrix
#' @author Chris Wallace, Mary Fortune
#' @export
ldak.weights <- function(XX,subset=NULL) {
    shell.path <- system.file("bash",package="abcGWAS")
    shell.exec <- file.path(shell.path,"run-ldak.sh")
    d <- tempdir()
    if(!is.null(subset)) {
        allsnps <- colnames(XX)
        XX <- XX[,subset]
    }
    write.plink(file.base=file.path(d,"ldak"), snps=XX,
            phenotype=rep(1,nrow(XX)),
            pedigree=1:nrow(XX),
            id=rep(1,nrow(XX)),
            mother=rep(0,nrow(XX)),
            father=rep(0,nrow(XX)),
            sex=rep(1,nrow(XX)),
            chromosome=rep(1,ncol(XX)))
    c1 <- paste(shell.exec,d)
    system(c1)
    weight_table <- read.table(file.path(d,"weightsALL"))
    snpswithp <- colnames(XX)
    weights <- weight_table$V1[unlist(lapply(snpswithp,function(x){which(weight_table$V5 == x)}))]
    names(weights) <- colnames(XX)
    if(!is.null(subset)) {
        dropped <- setdiff(allsnps,subset)
        weights.dropped <- structure(rep(0,length(dropped)),names=dropped)
        weights <- c(weights,weights.dropped)[allsnps]
    }
    weights
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Calculate distribution of wsumsq(|Zstar|,|Zo|) with Zstar sampled about Zo
##' @param Zo observed Z score vector
##' @param snp.data reference snp data
##' @param weights weights computed by ldak.weights
##' @param n specify number of samples to generate distance distribution (optional)
##' @return vector of distances
##' @author Chris Wallace
##' @export
Zstar.dist <- function(Zo,snp.data,weights,n=100000) {
    use <- which(weights>0)
    XX <- snp.data[,use]
    LD <- snpStats::ld(XX,XX,stat="R",symmetric=TRUE)
    LD<-as.matrix(make.positive.definite(LD))
    Zstar<-rmvnorm(n=n,mean=Zo[use],sigma=2*LD)
    wsumsqmat(abs(Zo[use]),t(abs(Zstar)),weights[use])
}
