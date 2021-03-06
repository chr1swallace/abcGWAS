##' @title computes the qth quartile of the distribution of distances about the Observed Z Score, Zo, at the true causal model. 
##' @export
##' @param Zo The Observed Z Score at the SNPs we shall analyse
##' @param LD	The LD matrix 
##' @param q The quartile we wish to compute epsilon at
##' @param weights The weights we get from LDAK for each SNP
##' @return The value of espilon, the ABC threshold
##' @author Mary Fortune 
compute_eps<-function(Zo,LD,q,weights){
	Zstar<-rmvnorm(n=10000,mean=Zo,sigma=LD)
	dist<-apply(Zstar,1,function(Z){wsumsq(abs(Zo),abs(Z),weights)})
	eps<-quantile(ecdf(dist),q)
	return(eps)
}

