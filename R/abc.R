##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title generate ABC samples for single causal variant models
##' @param Zo observed Z score
##' @param snp.data SnpMatrix of reference data
##' @param weights LDAK weights
##' @param freq phased reference haplotypes with Probability column
##' @param eps maximum difference between simulated and observed data
##'     for samples to pass
##' @param outdir directory within which output file will be created
##' @param nsam_model number of causal models to consider
##' @param nsam_gamma number of gamma to consider per model
##' @return filename within which samples that passed the threshold
##'     are written
##' @author Chris Wallace, Mary Fortune
##' @importFrom simGWAS make_GenoProbList
##' @export
run1CV <- function(Zo,snp.data,weights,freq,eps,outdir,nsam_model=100,nsam_gamma=1000) {
    snps <- colnames(snp.data)
    snps.use <- names(weights)[weights>0]
    if(!all(snps.use %in% snps))
        stop("snps with weights must be a subset of snp.data")
    if(!identical(c(snps,"Probability"),colnames(freq)))
        stop("all snps must be represented in freq")
    
    ## subset of SNPs at which sumsq will be evaluated
    weights.use<-weights[snps.use]
    snp.data.use <- snp.data[,snps.use]
    Zo <- Zo[snps.use]
    
    ## LD which will be needed for Zsim
    LD <- snpStats::ld(snp.data.use,snp.data.use,stat="R",symmetric=TRUE)
    LD <- as.matrix(make.positive.definite(LD))

    ## output
    outfile<-tempfile(tmpdir =outdir,fileext=".out")
    cat("CV\tgamma\tsumsq\n",file=outfile)
    
    totest<-sample(1:length(snps),nsam_model)
    numtested<-rep(0,length(snps))
    ## loop over CVs to test
    pb <- txtProgressBar(min = 0, max = nsam_model, style = 3)
    i <- 1
    for(which.CV.test in totest) {
        setTxtProgressBar(pb, i); i <- i+1
        CV.test<-snps[which.CV.test]
        GenoProbList<-make_GenoProbList(snps.use,CV.test,freq)
    ##    R2_true_CV<-LD_all[which.CV,which.CV.test]^2
        gamma1.test<-rgamma(nsam_gamma,0.8,1.2)*sample(c(-1,1),nsam_gamma,replace=TRUE)
        Zsim <- sapply(gamma1.test, function(gamma1) {
            Ze <- est_statistic(N0,N1,snps.use,CV.test,gamma1,freq,GenoProbList)            
            rmvnorm(n=1,mean=Ze,sigma=LD)
        })
        sumsq <- wsumsqmat(abs(Zo),abs(Zsim),weights.use)
        pass <- which(sumsq<eps)
        if(length(pass)) {
            ret <- data.frame(CV.test,gamma1.test[pass],sumsq[pass])
            write.table(ret,file=outfile,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
        }
    }
    close(pb)
    for(i in totest)
        numtested[i] <- numtested[i] + nsam_gamma
    save(numtested,file=sub(".out$",".RData",outfile))
    return(outfile)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title generate ABC samples for single causal variant models
##' @param Zo observed Z score
##' @param snp.data SnpMatrix of reference data
##' @param weights LDAK weights
##' @param freq phased reference haplotypes with Probability column
##' @param eps maximum difference between simulated and observed data
##'     for samples to pass
##' @param outdir directory within which output file will be created
##' @param nsam_model number of causal models to consider
##' @param nsam_gamma number of gamma to consider per model
##' @return filename within which samples that passed the threshold
##'     are written
##' @author Chris Wallace, Mary Fortune
##' @export
abcfillgaps <- function(Zo,snp.data,weights,freq,eps,outdir,nsam_model=100,nsam_gamma=1000,ncv=1,
                        uniqify=TRUE,method=c("random","gapfill")) {
    snps <- colnames(snp.data)
    snps.use <- names(weights)[weights>0]
    if(!all(snps.use %in% snps))
        stop("snps with weights must be a subset of snp.data")
    if(!identical(c(snps,"Probability"),colnames(freq)))
        stop("all snps must be represented in freq")
    
    ## subset of SNPs at which sumsq will be evaluated
    weights.use<-weights[snps.use]
    snp.data.use <- snp.data[,snps.use]
    Zo <- Zo[snps.use]
    
    ## LD which will be needed for Zsim
    LD <- snpStats::ld(snp.data.use,snp.data.use,stat="R",symmetric=TRUE)
    LD <- as.matrix(make.positive.definite(LD))
    
    ## output
    outfile<-tempfile(tmpdir =outdir,fileext=".out")
    cat("CV\tgamma\tsumsq\n",file=outfile)

    method <- match.arg(method)
    if(method=="random") {
        if(ncv==1) {
            totest<-sample(1:length(snps),nsam_model)
        } else {
            totest <- replicate(nsam_model,sort(sample(1:length(snps),ncv)),simplify=FALSE)
            if(uniqify) {
                dups <- which(duplicated(sapply(totest,paste,collapse="-")))
                if(length(dups))
                    totest <- totest[-dups]
            }
        }
    } else {
        sfile <- summarise.coverage(outdir)
        (load(sfile))
        if(ncv==1) {
            totest <- which(nt1==0)
        } else {
            wh <- which(nt2==0,arr.ind=TRUE)
            wh <- wh[ wh[,1] < wh[,2], ]
            wh <- t(wh)
            totest <- lapply(1:ncol(wh),function(i) wh[,i])
        }        
        if(!length(totest))
            stop("no more samples required for ncv =",ncv)
    }
    nt1<-rep(0,length(snps))
    nt2<-matrix(0,length(snps),length(snps))
    ## loop over CVs to test
    pb <- txtProgressBar(min = 0, max = length(totest), style = 3)
    for(i in seq_along(totest)) {
        setTxtProgressBar(pb, i); 
        which.CV.test <- totest[[i]]
        CV.test<-snps[which.CV.test]
        cv.str <- paste(CV.test,collapse=" ")
        GenoProbList<-make_GenoProbList(snps.use,CV.test,freq)
        ##    R2_true_CV<-LD_all[which.CV,which.CV.test]^2
        gamma1.test<-rgamma(nsam_gamma*ncv,0.8,1.2)*sample(c(-1,1),nsam_gamma*ncv,replace=TRUE)
        if(ncv>1) {
            gamma1.test<-split(gamma1.test,rep(1:nsam_gamma,ncv))
            gamma1.str <- sapply(gamma1.test,paste,collapse=" ")
        } else {
            gamma1.str <- gamma1.test
        }
        Zsim <- sapply(seq_along(gamma1.test), function(i) {
            Ze <- est_statistic(N0,N1,snps.use,CV.test,gamma1.test[[i]],freq,GenoProbList)            
            rmvnorm(n=1,mean=Ze,sigma=LD)
        })
        sumsq <- wsumsqmat(abs(Zo),abs(Zsim),weights.use)
        pass <- which(sumsq<eps)
        if(length(pass)) {
            ret <- data.frame(cv.str,gamma1.str[pass],sumsq[pass])
            write.table(ret,file=outfile,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
        }
    }
    close(pb)    
    if(ncv==1) {
        for(i in totest)
            nt1[i] <- nt1[i] + nsam_gamma
    } else {
        for(i in seq_along(totest))
            nt2[totest[[i]][1],totest[[i]][2]] <- nt2[totest[[i]][1],totest[[i]][2]] + 1
    }
    save(nt1,nt2,file=sub(".out$",".RData",outfile))
    return(outfile)
}
##' @title summarise coverage of an ABC run to date
##' @param d directory containing abc samples
##' @param force default FALSE, in which case summary.RData will be regenerated only if one of the abc sample files is newer. If TRUE, regerate anyway
##' @return filename of summary file
##' @export
##' @author Chris Wallace
summarise.coverage <- function(d,force=FALSE) {
    files <- list.files(d,pattern="file.*.RData",full=TRUE)
    ofile <- file.path(d,"summary.RData")
    if(!force & file.exists(ofile)) {
        i.mtime <- min(file.mtime(files))
        o.mtime <- file.mtime(ofile)
        if(o.mtime>i.mtime) {
            message("summary file already up to date, use force=TRUE to regenerate")
            return(ofile)
        }
    } 
    NT <- model.coverage(files)
    nt1 <- NT$nt1
    nt2 <- NT$nt2
    counts2 <- nt2[upper.tri(nt2)]
    cat(min(nt1),min(counts2),file=sub(".RData",".txt",ofile))
    save(nt1,nt2,file=ofile)
    return(ofile)
}
model.coverage <- function(files) {
        load(files[1])
        NT1 <- nt1
        NT2 <- nt2
        if(length(files)>1) {
            for(i in 2:length(files)) {
                (load(files[i]))
                NT1 <- NT1 + nt1
                NT2 <- NT2 + nt2
            }
        }
        message("single CV models")
        summary(NT1)
        message("two CV model visits")
        summary(counts2 <- NT2[upper.tri(NT2)])
        return(list(nt1=NT1,nt2=NT2))
}
    
