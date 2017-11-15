
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title generate ABC samples for causal variant models
##' @param Zo observed Z score
##' @param snp.data SnpMatrix of reference data
##' @param weights LDAK weights
##' @param freq phased reference haplotypes with Probability column
##' @param eps maximum difference between simulated and observed data
##'     for samples to pass
##' @param outdir directory within which output file will be created
##' @param N0 number of controls
##' @param N1 number of cases
##' @param nsam_model number of causal models to consider
##' @param nsam_gamma number of gamma to consider per model
##' @param ncv number of causal variants
##' @param uniqify if TRUE, don't sample the same model more than nsam_gamma times
##' @param method select models at random, or to fill gaps in already sampled
##' @return filename within which samples that passed the threshold
##'     are written
##' @author Chris Wallace, Mary Fortune
##' @export
abcsample <- function(Zo,snp.data,weights,freq,eps,outdir,N0,N1,nsam_model=100,nsam_gamma=1000,ncv=1,
                        uniqify=TRUE,method=c("random","gapfill")) {
    snps <- setdiff(colnames(freq),"Probability")
    snps.use <- names(weights)[weights>0]
    if(!all(snps.use %in% snps))
        stop("snps with non-zero weights must be a subset of freq data")
    ## if(!identical(c(snps,"Probability"),colnames(freq)))
    if(!all(snps.use %in% colnames(snp.data)))
        stop("all snps with non-zero weights must be represented in snp.data")
    if(is.null(names(Zo)))
        stop("Zo is unnamed, must have names that are a subset of colnames(snp.data)")
    if(!all(snps.use %in% names(Zo)))
        stop("all snps with non-zero weights must have be represented in Zo")

    message("SNPs in reference data: ",length(snps))
    message("SNPs in observed Z score: ",length(Zo))
    message("SNPs with non-zero weights: ",length(snps.use))
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
    sfile <- summarise_coverage(outdir,filename.only=TRUE)
    if(method=="random"|| !file.exists(sfile)) {
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
        sfile <- summarise_coverage(outdir,force=TRUE)
        (load(sfile))
        if(ncv==1) {
            mincov <- min(nt1)
            totest <- which(nt1==mincov)
        } else {
            mincov <- min(nt2[upper.tri(nt2)])
            wh <- which(nt2==mincov,arr.ind=TRUE)
            wh <- wh[ wh[,1] < wh[,2], ]
            if(!nrow(wh))
                stop("no more samples required for ncv =",ncv)
            wh <- t(wh)
            totest <- lapply(1:ncol(wh),function(i) wh[,i])
        }        
        if(!length(totest))
            stop("no more samples required for ncv =",ncv)
        if(length(totest)>nsam_model)
            totest <- totest[sample(1:length(totest),nsam_model)]
    }
    nt1<-rep(0,length(snps))
    nt2<-matrix(0,length(snps),length(snps))
    ## loop over CVs to test
    pb <- txtProgressBar(min = 0, max = length(totest), style = 3)
    for(i in seq_along(totest)) {
        setTxtProgressBar(pb, i); 
        which.CV.test <- totest[[i]]
        ## print(totest[[i]])
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
            nt2[totest[[i]][1],totest[[i]][2]] <- nt2[totest[[i]][1],totest[[i]][2]] + nsam_gamma
    }
    save(nt1,nt2,file=sub(".out$",".RData",outfile))
    return(outfile)
}


##' @title summarise_coverage of an ABC run to date
##' @param d directory containing abc samples
##' @param force default FALSE, in which case summary.RData will be regenerated only if one of the abc sample files is newer. If TRUE, regerate anyway
##' @return filename of summary file
##' @author Chris Wallace
##' @export
summarise_coverage <- function(d,filename.only=FALSE,force=FALSE) {
    ofile <- file.path(d,"summary.RData")
    if(filename.only)
        return(ofile)
    files <- list.files(d,pattern="file.*.RData",full=TRUE)
    if(!length(files)) 
        stop("no sample output files found in ",d)
    if(!force & file.exists(ofile)) {
        i.mtime <- max(file.mtime(files))
        o.mtime <- file.mtime(ofile)
        if(o.mtime>i.mtime) {
            message("summary file already up to date, use force=TRUE to regenerate")
            return(ofile)
        }
    }
    NT <- model.coverage(files)
    nt1 <- NT$nt1
    nt2 <- NT$nt2
    if(all(nt2==0))
        nt2 <- NULL
    if(!is.null(nt2)) {
        counts2 <- nt2[upper.tri(nt2)]
    } else {
        counts2 <- 0
    }
    cat(min(nt1),min(counts2),file=sub(".RData",".txt",ofile))
    save(nt1,nt2,file=ofile)
    ## message("1 CV models: ",nt1)
    ## message("2 CV models: ",nt2)
    invisible(ofile)
}
model.coverage <- function(files) {
    obj <- load(files[1])
    NT1 <- nt1
    if("nt2" %in% obj) {        
        NT2 <- nt2
    } else {  ## for some old files, nt2 not saved if only 1 CV models run
        NT2 <- NULL
    }
    if(length(files)>1) {
        for(i in 2:length(files)) {
            (load(files[i]))
            NT1 <- NT1 + nt1
            if("nt2" %in% obj) {
                NT2 <- NT2 + nt2
            }
        }
    }
    ## message("single CV models")
    ## summary(NT1)
    ## message("two CV model visits")
    ## summary(counts2 <- NT2[upper.tri(NT2)])
    return(list(nt1=NT1,nt2=NT2))
}
##' @title print a summary of model coverage to date
##' @param d directory
##' @return list of nt1 and nt2, as loaded from summary file (invisibly returned)
##' @author Chris Wallace
##' @export
show.summary <- function(d) {
    (load(summarise_coverage(d,filename.only = FALSE)))
    message("single cv models:")
    tt <- table(nt1)
    for(n in names(tt))
        message("\t",n," samples:\t",tt[n])
    message("two cv models:")
    if(is.null(nt2)) {
        message("none found")
    } else {
        tt <- table(nt2)
        for(n in names(tt))
            message("\t",n," samples:\t",tt[n])
    }
    invisible(list(nt1=nt1,nt2=nt2))
}
    

################################################################################

## old stuff to be deleted

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

