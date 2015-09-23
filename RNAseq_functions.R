######################################################
### Cas Retel
### Functions used during RNAseq research
### VUmc & Leiden University, 02-2014 until 07-2015
######################################################

#############################
### Main predictive functions, as discussed in thesis
predictPower <- function(nsample, b0, b1, size, q0, delta=0, alpha=0.1, 
                         fixedseed=FALSE){
  # Computes power and effect abundance for the detection of differentially 
  # expressed features in a simulated dataset. Requires 
  # (possibly zero-inflated) negative binomial coefficients as input, 
  # generates a dataset and performs exact tests of differential expression 
  # based on (Robinson, 2007). Assumes an equalsized two-group comparison.  
  if(0 %in% sapply(list(nsample, b0, b1, size, q0, delta, alpha), 
                   function(x) length(x))){
    stop("Not all input arguments are specified")
  }
  if(any(is.na(c(nsample, b0, b1, size, q0)))) stop("Input contains NA values")
  if(!all.equal(length(b0), length(b1), length(size)) | 
       (!is.null(q0) & length(b0)!=length(q0))){
    stop("Coefficient vectors are not of the same length")
  }
  require(edgeR)
  p <- length(b0)
  
  # round to nearest smaller even number: 
  while(nsample%%2) nsample <- ceiling(nsample - 1)  
  if(is.null(q0)) q0 <- numeric(p)  # when likelihood is nb
  
  ### Generate count matrix:  
  if(fixedseed){
    if(!exists(".Random.seed", mode="integer", envir=globalenv())) sample(NA)
    oldSeed <- .Random.seed
    set.seed(2909)
  }
  sd <- matrix(0, p, nsample)  # matrix of counts
  sd[, 1:(nsample/2)] <- t(sapply(1:p, function(i) 
    rnbinom(n=nsample/2, size=exp(size[i]), mu=exp(b0[i]))))
  sd[, ((nsample/2)+1):nsample] <- t(sapply(1:p, function(i) 
    rnbinom(n=nsample/2, size=exp(size[i]), mu=exp(b0[i]+b1[i]))))
  for(i in 1:p) sd[i, sample(1:nsample, size=round(nsample*q0[i]))] <- 0
  set.seed(oldSeed)
  
  ### Exact testing:
  res <- DGEList(sd, norm.factors=calcNormFactors(sd, method="TMM"), 
                 group=factor(rep(1:2, each=ncol(sd)/2)))
  res <- estimateCommonDisp(res)
  res <- estimateTagwiseDisp(res)
  
  res <- exactTest(res)
  res <- topTags(res, n=nrow(sd), sort.by="none")$table
  
  ### Output
  out <- data.frame(ActualPos=sum(b1!=0), Rejections=sum(res$FDR <= alpha), 
                    TruePos=sum(which(b1!=0) %in% which(res$FDR <= alpha)))
  #   out <- data.frame(DE.feats=sum(b1!=0), detected.feats=sum(res$FDR <= alpha), 
  #                     matching.feats=sum(is.element(which(res$FDR <= alpha), 
  #                                                   which(b1!=0))))
  if(delta!=0){
    out$ActualPos.delta=sum(abs(b1)>=delta)
    out$TruePos.delta=
      sum(which(res$FDR <= alpha) %in% which(abs(b1)>=delta))
  }
  rm(sd, res, p)
  return(out)
}


powerCurve <- function(fasobj, muprobj, nsample=c(10, 20, 40, 60, 80, 100), 
                       delta=0, alpha=0.1, niter=50, ncpus=2){
  # Computes power and effect abundance for the detection of 
  # differentially expressed features in a simulated dataset of 
  # negative binomial counts. Requires (possibly zero-inflated) 
  # negative binomial coefficients as input, generates a dataset and 
  # performs exact tests of differential expression based on (Robinson, 2007). 
  # Assumes an equalsized two-group comparison. 
  require(edgeR); require(INLA); require(snowfall)
  b1s <- function(p, obj, fixednumbernon0=TRUE){
    # draw p samples from updated prior distribution
    # fixednumbernon0 = should the number of nonzero effects be fixed
    #                   or determined stochastically?
    
    if(!any(is.element(c("p0est", "best"), names(obj)))){
      stop("input is not an UpdatePrior-object")
    }
    effect <- numeric(p)
    p0 <- effectindex <- numeric()
    
    p0 <- ifelse(is.element("p0est", names(obj)), obj$p0est, obj$best["p0"])
    
    if(fixednumbernon0){
      effectindex <- sample(1:p, size=round((1-p0)*p))
    }else{
      effectindex <- runif(p) > p0
    }
    
    # for NonPara-objects:
    if("p0est" %in% names(obj)){    # for NonParaUpdate-objects: 
      effect[effectindex] <- inla.rmarginal(n=sum(effectindex != 0), 
                                            obj$priornew)
    }else{    # for MixtureUpdate-objects:
      effect[effectindex] <- rnorm(sum(effectindex!=0), mean=obj$best["mu"], 
                                   sd=obj$best["stdev"])
      effect[effect!=0] <- effect[effect!=0] + 
        (-2 * (runif(sum(effectindex!=0)) > obj$best["pmin"]) * 
           effect[effect!=0])
      # because x + (-2*x) = -x
    }
    return(effect)
  }
  
  predictPower <- function(nsample, b0, b1, size, q0, delta=0, alpha=0.1, 
                           fixedseed=fixedseed){
    
    if(0 %in% unlist(lapply(list(nsample, b0, b1, size, q0, delta, alpha), 
                            function(x) length(x)))){
      stop("Not all input arguments are specified")
    }
    if(any(is.na(c(nsample, b0, b1, size, q0)))) stop("Input contains NA values")
    if(!all.equal(length(b0), length(b1), length(size)) | 
         (!is.null(q0) & length(b0)!=length(q0))){
      stop("Coefficient vectors are not of the same length")
    }
    require(edgeR)
    p <- length(b0)
    
    # round to nearest smaller even number: 
    while(nsample%%2) nsample <- ceiling(nsample - 1)  
    if(is.null(q0)) q0 <- numeric(p)  # when likelihood is nb
    
    ### Generate count matrix:  
    if(fixedseed){
      if(!exists(".Random.seed", mode="integer", envir=globalenv())) sample(NA)
      oldSeed <- .Random.seed
      set.seed(2909)
    }
    sd <- matrix(0, p, nsample)  # matrix of counts
    sd[, 1:(nsample/2)] <- t(sapply(1:p, function(i) 
      rnbinom(n=nsample/2, size=exp(size[i]), mu=exp(b0[i]))))
    sd[, ((nsample/2)+1):nsample] <- t(sapply(1:p, function(i) 
      rnbinom(n=nsample/2, size=exp(size[i]), mu=exp(b0[i]+b1[i]))))
    for(i in 1:p) sd[i, sample(1:nsample, size=round(nsample*q0[i]))] <- 0
    set.seed(oldSeed)
    
    ### Exact testing:
    res <- DGEList(sd, norm.factors=calcNormFactors(sd, method="TMM"), 
                   group=factor(rep(1:2, each=ncol(sd)/2)))
    res <- estimateCommonDisp(res)
    res <- estimateTagwiseDisp(res)
    
    res <- exactTest(res)
    res <- topTags(res, n=nrow(sd), sort.by="none")$table
    
    ### Output
    out <- data.frame(ActualPos=sum(b1!=0), Rejections=sum(res$FDR <= alpha), 
                      TruePos=sum(which(b1!=0) %in% which(res$FDR <= alpha)))
    #   out <- data.frame(DE.feats=sum(b1!=0), detected.feats=sum(res$FDR <= alpha), 
    #                     matching.feats=sum(is.element(which(res$FDR <= alpha), 
    #                                                   which(b1!=0))))
    if(delta!=0){
      out$ActualPos.delta=sum(abs(b1)>=delta)
      out$Truepos.delta=
        sum(which(res$FDR <= alpha) %in% which(abs(b1)>=delta))
    }
    rm(sd, res, p)
    return(out)
  }
  
  p <- length(fasobj$res)
  b0 <- extB0s(fasobj)
  size <- extSizes(fasobj)
  q0 <- try(extCountp0s(fasobj))
  
  sfInit(cpus=ncpus, parallel=TRUE)
  sfLibrary(edgeR); sfLibrary(INLA)
  
  sfExport("predictPower", "b1s", "p","nsample", "niter", 
           "alpha", "b0", "size", "q0", "muprobj", "delta")
  
  abundance <- lapply(seq_along(nsample), function(i) matrix(NA, niter, 5))
  abundance <- lapply(nsample, function(obj)
    t(sfSapply(1:niter, function(i)
      predictPower(obj, b0, b1s(p, muprobj), size, q0, delta=delta, alpha=alpha))))
  
  sfStop()
  
  out <- t(sapply(abundance, function(x) 
    apply(x, 2, function(y) mean(as.numeric(y)))))
  out <- round(out)
  out <- cbind(nsample, out)
  
  colnames(out) <- c("nsample","ActualPos", "Rejections", "TruePos", 
                     "ActualPos.delta", "TruePos.delta")[1:ncol(out)]
  
  rm(abundance, p, b0, size, q0)
  if(b1thresh!=0) return(out) else return(out[, 1:4])
}



#############################
### ShrinkBayes convenience functions
extB0s <- function(...){
  # ExtB0s() returns estimates beta0 (intercept) coefficients 
  # from input FitAllShrink-objects
  # Returns na when when a $res-list is NULL
  fasall <- list(...)
  if(any(!sapply(fasall, function(obj) 
    ("res" %in% names(obj) | "priors" %in% names(obj))))){
    stop("Input are not all FitAllShrink-objects")
  }
#   if(any(!sapply(seq_along(fasall), function(i) 
#     ("res" %in% names(fasall[[i]]) | "priors" %in% names(fasall[[i]]))))){
#     stop("Input are not all FitAllShrink-objects")
#   }
  as.numeric(sapply(fasall, function(obj) 
    sapply(1:length(obj$res), function(i) as.numeric(
      obj$res[[i]]$summary.fixed[1, 1])))) 
}

extB1s <- function(...){ 
  # ExtB1s() returns estimated beta1 (group effect) coefficients 
  # from FitAllShrink-objects
  # Returns na when when a $res-list is NULL
  fasall <- list(...)
  if(any(!sapply(fasall, function(obj) 
    ("res" %in% names(obj) | "priors" %in% names(obj))))){
    stop("Input are not all FitAllShrink-objects")
  }
  as.numeric(sapply(fasall, function(obj) 
    sapply(1:length(obj$res), function(i) as.numeric(
      obj$res[[i]]$summary.fixed[2, 1])))) 
}
extSizes <- function(..., logscale=T){ 
  # ExtSizes() returns estimated log(size) coefficients 
  # from FitAllShrink-objects. 
  # set logscale = F to return the exponent of this estimate
  # Returns na when when a $res-list is NULL
  fasall <- list(...)
  if(any(!sapply(fasall, function(obj) 
    ("res" %in% names(obj) | "priors" %in% names(obj))))){
    stop("Input are not all FitAllShrink-objects")
  }
  out <- as.numeric(sapply(fasall, function(obj) 
    sapply(1:length(obj$res), function(i) as.numeric(
      obj$res[[i]]$internal.summary.hyper[1, 1]))))
  if(!logscale) out <- exp(out)
  return(out)
}
extCountp0s <- function(...){ 
  # extCountp0s() returns count null probabilities (q0)
  # from FitAllShrink-objects
  # Returns na when when a $res-list is NULL
  fasall <- list(...)
  if(any(!sapply(fasall, function(obj) 
    ("res" %in% names(obj) | "priors" %in% names(obj))))){
    stop("Input are not all FitAllShrink-objects")
  }
  as.numeric(sapply(fasall, function(obj) 
    sapply(1:length(obj$res), function(i) as.numeric(
      obj$res[[i]]$summary.hyperpar[2, 1]))))
}
extPriorEffp0s <- function(..., digits=NA){
  # Extracts effect null probabilities from UpdatePrior-objects 
  # (either Mixture or NonPara). Returns a numeric vector. 
  # digits is an optional parameter, specifying number of digits to round to
  objall <- list(...)
  
  if(any(!sapply(objall, function(obj) 
    ("p0est" %in% names(obj) | "best" %in% names(obj))))){
      stop("input is not an UpdatePrior-object")
  }
    
  out <- numeric(length(objall))
  out <- sapply(seq_along(objall), function(i) 
    if(any(names(objall[[i]])=="best")){  # Mixture prior
      out <- objall[[i]]$best["p0"] 
    } else out <- objall[[i]]$p0est)   # Nonparametric prior
  if(!is.na(digits)) out <- round(out, digits=digits)
  as.numeric(out)
}

extEffp0s <- function(obj){
  # Extracts effect null probabilities from an UpdatePosterior-object
  # (either Mixture or NonPara). Returns a numeric vector. 
  if(!identical(names(obj[[1]]), c("postbetanon0", "postbeta0", "loglik"))){
    stop("input is not an UpdatePosterior-object")
  }
  
  as.numeric(sapply(obj, function(x) x$postbeta0))
}

sampleB1s <- function(obj, expectation=FALSE, numbernon0=NA){
  # sampleB1s() returns a draw from posterior distributions in the
  # input UpdatePosterior-object (either Mixture or NonPara). 
  # expectation = Boolean. Return expected value instead of draw?
  # numbernon0 = Numeric. When specified, those numbernon0 distributions
  #              with smallest posterior p0 are used to draw from 
  if(is.na(numbernon0)){   # zero or draw determined stochastically
    effectind <- (runif(length(obj)) > 
                    sapply(obj, function(x) x$postbeta0))
  }else{   # take smallest p0's
    effectind <- is.element(
      seq_along(obj), order(unlist(
        lapply(obj, function(x) x$postbeta0)))[1:numbernon0])
  }
  
  effect <- numeric(length(obj))
  if(!expectation){   # random draw
    effect[effectind] <- sapply(which(effectind), function(i) 
      inla.rmarginal(1, obj[[i]]$postbetanon0[[1]]))
  }else{   # expected value
    effect[effectind] <- sapply(which(effectind), function(i) 
      inla.emarginal(function(x) x, obj[[i]]$postbetanon0[[1]]))
  }
  return(effect)
}

expB1s <- function(obj, ignorep0=FALSE, stochasticp0=FALSE){
  # expB1s() returns expected values of (posterior distribution times the
  # corresponding null probability), for every feature, from an 
  # UpdatePosterior-object. 
  # Set ignorep0 = TRUE to return expectation of nonzero part only
  # Set stochasticp0 = TRUE to compare effect p0 to a random draw 
  # (with fixed seed) from U(0, 1) distr., and use this to determine 
  # whether 0 or E(nonzero dist) is returned. 
  if(!is.list(obj) | 
       !identical(names(obj[[1]]), c("postbetanon0", "postbeta0", "loglik"))){
    stop("input is not an UpdatePosterior-object")
  }
  p <- length(obj)  
  effect <- numeric(p)
  
  # return 0 for posteriors with zero probability of exactly 1: 
  # for these posteriors, marginal INLA functions will result in errors
  nonzeroind <- which(sapply(obj, function(x) (x$postbeta0!=1)))
  
  if(ignorep0){  # return E(nonzero part)
    effect[nonzeroind] <- sapply(nonzeroind, function(i) 
      inla.emarginal(function(x) x, obj[[i]]$postbetanon0[[1]]))
  }else{
    if(!stochasticp0){  # effectp0*E(nonzero part)
      effect[nonzeroind] <- sapply(nonzeroind, function(i) 
        (1-obj[[i]]$postbeta0) * 
          inla.emarginal(function(x) x, obj[[i]]$postbetanon0[[1]]))
    }else{  # return (pseudo)random sample from E(nonzero)'s
      #if(!exists(".Random.seed", mode="integer", envir=globalenv())) sample(NA)
      #oldSeed <- .Random.seed
      #set.seed(1605)  
      
      nonzeroind <- runif(p) > sapply(1:p, function(i) obj[[i]]$postbeta0)
      effect[nonzeroind] <- sapply(which(nonzeroind), function(i) 
        inla.emarginal(function(x) x, obj[[i]]$postbetanon0[[1]]))
      
      #assign(".Random.seed", oldSeed, envir=globalenv())
    }
  }
  return(effect)
}

samplePriorSizes <- function(p, obj, logscale=T){
  # draws n samples from prior size distribution, as stored in ShrinkSeq-object
  if(logscale){
    return(rnorm(p, mean=obj$pmlist$mudisp, sd=1/sqrt(obj$pmlist$precdisp)))
  }else{
    return(exp(rnorm(p, mean=obj$pmlist$mudisp, 
                     sd=1/sqrt(obj$pmlist$precdisp))))
  }
}

samplePriorB1s <- function(p, obj, fixednumbernon0=TRUE){
  # draw p samples from updated prior distribution
  # fixednumbernon0 = should the number of nonzero effects be fixed
  #                   or determined stochastically?
  
  if(!any(is.element(c("p0est", "best"), names(obj)))){
    stop("input is not an UpdatePrior-object")
  }
  effect <- numeric(p)
  p0 <- effectindex <- numeric()
  
  p0 <- ifelse(is.element("p0est", names(obj)), obj$p0est, obj$best["p0"])
  
  if(fixednumbernon0){
    effectindex <- sample(1:p, size=round((1-p0)*p))
  }else{
    effectindex <- runif(p) > p0
  }
  
  # for NonPara-objects:
  if("p0est" %in% names(obj)){    # for NonParaUpdate-objects: 
    effect[effectindex] <- inla.rmarginal(n=sum(effectindex != 0), 
                                          obj$priornew)
  }else{    # for MixtureUpdate-objects:
    effect[effectindex] <- rnorm(sum(effectindex!=0), mean=obj$best["mu"], 
                                 sd=obj$best["stdev"])
    effect[effect!=0] <- effect[effect!=0] + 
      (-2 * (runif(sum(effectindex!=0)) > obj$best["pmin"]) * 
         effect[effect!=0])
    # because x + (-2*x) = -x
  }
  return(effect)
}


#############################
### Plotting of ShrinkBayes output
plotMixturePrior <- function(mixtureupdate, scaletodens=TRUE, plotp0=TRUE, 
                             colour=1, ltype=1, lwdth=1, xlim=c(-5, 5), 
                             ylim=NA, main=NULL){
  # scaletodens = logical indicating whether the nonzero part of the 
  # distribution (that integrates to 1) should be scaled, so that
  # [ p0*diracdelta(0) + (nonzero part) ] integrates to 1 
  # plotp0 = logical indicating whether the point mass on zero is plotted
  p0 <- pneg <- ppos <- mu <- sd <- numeric(1)
  x <- y <- numeric(1000)
  
  p0 <- mixtureupdate$best["p0"]
  pneg <- mixtureupdate$best["pmin"]
  ppos <- (1-mixtureupdate$best["pmin"])
  if(scaletodens){
    pneg <- (1-p0)*pneg
    ppos <- (1-p0)*ppos
  }
  
  mu <- mixtureupdate$best["mu"]
  sd <- mixtureupdate$best["stdev"]
  
  x <- seq(xlim[1], xlim[2], length.out=1000)
  y <- pneg*dnorm(x, -mu, sd) + ppos*dnorm(x, mu, sd)
  suppressWarnings(if(is.na(ylim)) ylim <- c(0, 1.1*max(y, p0)))
  plot(x, y, type='l', col=colour, lty=ltype, lwd=lwdth, ylim=ylim, 
       main=main, ylab="Pr(x)")
  if(plotp0) points(0, p0, pch=(15+ltype), col=colour)
}


linesMixturePrior <- function(mixtureupdate, scaletodens=TRUE, plotp0=TRUE, 
                              colour=1, ltype=1, lwdth=2, xlim=c(-5, 5)){
  # linesMixturePrior is meant to be used in conjunction with plotMixturePrior()
  p0 <- pneg <- ppos <- mu <- sd <- numeric(1)
  x <- y <- numeric(1000)
  
  p0 <- mixtureupdate$best["p0"]
  pneg <- mixtureupdate$best["pmin"]
  ppos <- (1-mixtureupdate$best["pmin"])
  if(scaletodens){
    pneg <- (1-p0)*pneg
    ppos <- (1-p0)*ppos
  }
  
  mu <- mixtureupdate$best["mu"]
  sd <- mixtureupdate$best["stdev"]
  
  x <- seq(xlim[1], xlim[2], length.out=1000)
  y <- pneg*dnorm(x, -mu, sd) + ppos*dnorm(x, mu, sd)
  lines(x, y, type='l', lty=ltype, col=colour, lwd=lwdth)
  if(plotp0) points(0, p0, pch=(15+ltype), col=colour)
}


plotNonParaPrior <- function(nonparaupdate, scaletodens=TRUE, plotp0=TRUE, 
                             colour=1, ltype=1, lwdth=1, xlim=c(NA, NA), 
                             ylim=NA, main=NULL){
  # scaletodens = logical indicating whether the nonzero part of the 
  # distribution (that integrates to 1) should be scaled, so that
  # [ p0*diracdelta(0) + (nonzero part) ] integrates to 1 
  # plotp0 = logical indicating whether the point mass on zero is plotted
  p0 <- x <- y <- numeric()
  
  p0 <- nonparaupdate$p0est
  x <- nonparaupdate$priornew[, "nx"]
  y <- nonparaupdate$priornew[, "ny"]
  if(scaletodens) y <- (1-p0)*y
  
  xmax <- ifelse(is.na(xlim[2]), max(x), xlim[2])
  xmin <- ifelse(is.na(xlim[1]), min(x), xlim[1])
  
  suppressWarnings(if(is.na(ylim)) ylim <- c(0, 1.1*max(y, p0)))
  plot(x, y, type='l', col=colour, lty=ltype, lwd=lwdth, xlim=c(xmin, xmax), 
       ylim=ylim, main=main, ylab="Pr(x)")
  if(plotp0) points(0, p0, pch=(15+ltype), col=colour)
}


linesNonParaPrior <- function(nonparaupdate, scaletodens=TRUE, plotp0=TRUE, 
                              colour=1, ltype=1, lwdth=2){
  # linesNonParaPrior is meant to be used in conjunction with plotNonParaPrior()
  p0 <- x <- y <- numeric()
  
  p0 <- nonparaupdate$p0est
  x <- nonparaupdate$priornew[, "nx"]
  y <- nonparaupdate$priornew[, "ny"]
  if(scaletodens) y <- (1-p0)*y
  
  lines(x, y, type='l', col=colour, lty=ltype, lwd=lwdth)
  if(plotp0) points(0, p0, pch=(15+ltype), col=colour)
}


plotFASPrior <- function(fitallshrink, xlimit=NA, ylimit=NA, colour=1, ltype=1){
  mu <- fitallshrink$priors$mufixed 
  sd <- sqrt(1/fitallshrink$priors$precfixed)
  suppressWarnings(if(is.na(xlimit)) xlimit <- c(mu-(4*sd), mu+(4*sd)))
  x <- seq(xlimit[1], xlimit[2], length.out=1000)
  y <- dnorm(x, mu, sd)
  suppressWarnings(if(is.na(ylimit)) ylimit <- c(0, 1.1*max(y)))
  plot(x, y, type='l', col=colour, lty=ltype, ylim=ylimit)
}


linesFASPrior <- function(fitallshrink, colour=1, ltype=1, xlimit=NA){
  # linesFASPrior is meant to be used in conjunction with plotFASprior()
  mu <- fitallshrink$priors$mufixed 
  sd <- sqrt(1/fitallshrink$priors$precfixed)
  suppressWarnings(if(is.na(xlimit)) xlimit <- c(mu-(6*sd), mu+(6*sd)))
  x <- seq(xlimit[1], xlimit[2], length.out=1000)
  y <- dnorm(x, mu, sd)
  lines(x, y, type='l', col=colour, lty=ltype, ylim=ylim)
}


# auc (= area under curve) is used to normalize the marginals obtained 
# from NonParaUpdatePosterior(), since these are not densities
auc <- function(dens){
  x <- dens[, 1]
  y <- dens[, 2]
  dx <- sapply(1:(length(x)-1), function(i) x[i+1]-x[i])
  sum(y[-1]*dx)
}

plotPosterior <- function(upobj, index=1, xlim=c(NA, NA), ymax=1.1){
  # index = numeric vector indicating which marginals to draw
  p <- length(index)
  xmax.data <- max(sapply(1:p, function(i) 
    max(upobj[[index[i]]]$postbetanon0[[1]][, "x"])))
  xmin.data <- min(sapply(1:p, function(i) 
    min(upobj[[index[i]]]$postbetanon0[[1]][, "x"])))
  xmax <- ifelse(is.na(xlim[2]), xmax.data, xlim[2])
  xmin <- ifelse(is.na(xlim[1]), xmin.data, xlim[1])
  with(upobj[[index[1]]], plot(
    x=postbetanon0[[1]][, "x"], 
    y=((1-postbeta0)/auc(postbetanon0[[1]]))*postbetanon0[[1]][, "y"], 
    type='l', xlim=c(xmin, xmax), ylim=c(0, ymax), xlab="beta", ylab="pr(beta)"
  ))
  with(upobj[[index[1]]], points(0, postbeta0, pch=16))
  if(p > 1){
    cols <- rainbow(p-1)
    for(i in 2:p){
      with(upobj[[index[i]]], lines(
        x=postbetanon0[[1]][, "x"], 
        y=((1-postbeta0)/auc(postbetanon0[[1]]))*postbetanon0[[1]][, "y"], 
        col=cols[i-1]
      ))
      with(upobj[[index[i]]], points(0, postbeta0, pch=16, col=cols[i-1])) 
    }
  }
}

linesPosterior <- function(upobj, index=1, lty=2, col=NA){
  # linesPosterior() is meant to be used in conjunction with plotPosterior();
  # same index will have same color, but different linetype. 
  # convenient to compare marginals calculated under different constraints
  p <- length(index)
  if(is.na(col)) col <- 1
  with(upobj[[index[1]]], lines(
    x=postbetanon0[[1]][, "x"], 
    y=((1-postbeta0)/auc(postbetanon0[[1]]))*postbetanon0[[1]][, "y"], 
    col=col, lty=lty
  ))
  with(upobj[[index[1]]], points(0, postbeta0, pch=16, col=col))
  if(p > 1){
    if(is.na(col)) cols <- rainbow(p-1) else cols=col
    for(i in 2:p){
      with(upobj[[index[i]]], lines(
        x=postbetanon0[[1]][, "x"], 
        y=((1-postbeta0)/auc(postbetanon0[[1]]))*postbetanon0[[1]][, "y"], 
        col=cols[i-1], lty=lty
      ))
      with(upobj[[index[i]]], points(0, postbeta0, pch=16, col=cols[i-1])) 
    }
  }  
}


#############################
### Data manipulation
edgeRnormalization <- function(counts){
  # edgeRnormalization normalizes a table of counts according to TMM-method, 
  # implemented in as is usually done when creating a DGEList-object
  require(edgeR)
  rellibsize <- colSums(counts)/exp(mean(log(colSums(counts))))
  nf = calcNormFactors(counts, method="TMM")*rellibsize
  normcounts = round(sweep(counts, 2, nf, "/"))
  normcounts
}

samp <- function(factor, npergroup=10, replace=FALSE){
  # samp() is used to take a random sample from a factor, with an equal
  # number of subjects from every factor level. output is an index vector. 
  # npergroup = size per factor level
  # replace = allow sampling with replacement?
  l <- levels(factor)
  ind <- matrix(NA, npergroup, length(l))
    
  if(any(sapply(seq_along(l), function(i) 
    length(which(factor==l[i])) < npergroup))){ 
    replace <- TRUE
    warning("factor has less entries than specified npergroup, 
            subjects will be recycled")
  }
  ind <- sapply(seq_along(l), function(i) 
    sample(which(factor==l[i]), size=npergroup, replace=replace))
  return(sort(as.numeric(ind)))
}

sampByFraction <- function(factor, n=20){
  # sampByFraction() is used to take a random sample from a factor, 
  # keeping the original ratio's of subjects per level intact
  # Output is an index vector
  # n = total sample size  
  l <- levels(factor)
  ntotal <- sapply(seq_along(l), function(i) sum(factor==l[i]))
  
  npergroup <- ntotal*(n/length(factor))
  
  # make sure the returned index is of length n:
  if(sum(round(npergroup)) < n) npergroup[which.max(npergroup%%1)] <- 
    npergroup[which.max(npergroup%%1)] + 1
  if(sum(round(npergroup)) > n) npergroup[which.min(npergroup%%1)] <- 
    npergroup[which.min(npergroup%%1)] - 1
  
  npergroup <- round(npergroup)
  
  ind <- sapply(seq_along(l), function(i) 
    sample(which(factor==l[i]), size=npergroup[i]))
  return(sort(as.numeric(unlist(ind))))
}

removeZeroes <- function(counts, propnon0=0.2, report=TRUE){
  # removeZeroes() removes all rows from a dataset that have less than 
  # propnon0 nonzero entries
  remove <- apply(counts, 1, function(row)
    sum(row!=0) < propnon0*ncol(counts))
  if(report) cat("removeZeroes():", sum(as.numeric(remove)), 
                 "out of", nrow(counts), "rows removed \n")
  return(counts[!remove, ])
}

createZeroes <- function(x, nzero){
  # assigns nzero elements of vector x the value 0
  # created for own convenience, slower than 
  #   x[sample(seq_along(x), size=nzero)] <- 0
  if(!is.numeric(nzero) | nzero < 0 | nzero > length(x)){
    stop("nzero must be positive and no larger than the length of x")
  }
  
  out <- numeric(length(x))
  index <- sample(seq_along(x), size=(length(x)-nzero))
  out[index] <- x[index]
  return(out)
} 

calcMode <- function(x, multiplemodes=T){
  # calcMode() returns the mode(s) of input vector x
  # set multiplemodes=F to return only the first element of the set of modes
  out <- as.numeric(names(table(x))[table(x) %in% max(table(x))])
  if(!multiplemodes) out <- out[1]
  cat("Largest occurence =", table(x)[table(x) %in% max(table(x))][1], "\n")
  return(out)
}


#############################
### processing data archives downloaded from TCGA data portal
TCGA.extFileNames <- function(TCGA_dir, mapping=c("gene", "exon", "spljxn"), 
                              technique="IlluminaHiSeq_RNAseq"){
  # extFileNames returns a list of all gene, exon or splice junction 
  # RNAseq filenames in a TCGA archive. 
  # TCGA_dir = path to unzipped .tgz file downloaded from TCGA data portal; 
  #   contains at least folders "METADATA" and "RNAseq"
  # mapping = type of expression signal 
  # technique = technical device used to acquire data
  #   (name of the directory above Level_3)
  #   ! These directories may not have consistent names, therefore they
  #   require creative input to resemble structure of the downloaded folder
  
  strng <- l3dir <- filenames <- character(0)
  strng <- paste(".", mapping[1], ".quant", sep="")
  l3dir <- paste(as.character(TCGA_dir), "/RNASeq/UNC__", 
                 technique, "/Level_3", sep="")
  filenames <- list.files(path=l3dir, pattern=strng)
  filenames
}

TCGA.filelistToBarcode <- function(filelist, form=c("UCSC", "complete")){
  # transforms names of tab-delimited files containing RNAseq data, 
  # as downloaded from TCGA data matrix, into TCGA subject barcodes 
  # set form = "UCSC" to abbreviate to the form "TCGA-xx-xxxx-xxA" 
  #              (used to match to UCSC clinical data files)
  
  n <- length(filelist)
  out <- numeric(n)
  out <- strsplit(filelist, "[.]")
  out <- sapply(1:n, function(i) out[[i]][2])
  if(form[1]=="UCSC") out <- substr(out, 1, 15)
  return(as.character(out))
}

TCGA.createCountMat <- function(
  TCGA_dir, filelist, unit=c("raw_counts", "median_length_normalized", "RPKM"), 
  technique="IlluminaHiSeq_RNAseq"
){
  # createCountMat() combines a set of tab delimited text files in a TCGA 
  # data archive, containing RNAseq data, into one data frame. 
  # Subjects are marked by their so-called TCGA-barcode (column names)
  # TCGA_dir = path to unzipped .tgz file downloaded from TCGA data portal; 
  #   contains at least folders "METADATA" and "RNAseq"  
  # filelist = a vector with names of text files, usually obtained from 
  #   TCGA.extFileNames() function
  # unit = unit to extract 
  
  n <- length(filelist)
  filenames <- character()
  out <- data.frame()
  
  filenames <- paste(as.character(TCGA_dir), "/RNASeq/UNC__", 
                     technique, "/Level_3/", filelist, sep="")
  out <- sapply(1:n, function(i) 
    as.numeric(read.delim(filenames[i], header=T)[, unit[1]]))
  rownames(out) <- read.delim(filenames[1], header=T)[, 1]
  colnames(out) <- TCGA.filelistToBarcode(filelist, "complete")
  return(out)
}

TCGA.checkFeatureNames <- function(
  TCGA_dir, filelist, technique="IlluminaHiSeq_RNAseq"
){
  # checkFeatureNames() returns TRUE when rownames of all files are identical, 
  # i.e. when every text file has information on the same feature 
  # Should be used prior to createCountMat, but is necessary only once 
  #    (hence a separate function)
  # TCGA_dir = path to unzipped .tgz file downloaded from TCGA data portal; 
  #   contains at least folders "METADATA" and "RNAseq"
  # filelist = a vector with names of text files, usually obtained from 
  #   TCGA.extFileNames() function
  n <- length(filelist)
  filenames <- character()
  featurenames <- data.frame()
  
  filenames <- paste(as.character(TCGA_dir), "/RNASeq/UNC__", 
                     technique, "/Level_3/", filelist, sep="")
  
  identical <- logical(0)
  featurenames <- sapply(1:n, function(i) 
    as.character(read.delim(filenames[i], header=T)[, 1]))
  identical <- sapply(2:n, function(i) 
    all.equal(featurenames[, 1], featurenames[, i]))
  return(!any(!identical))
}


#############################
### Obsolete power calculation functions; 
### used for exploration and analysis before calcPower() and powerCurve()
predictDEFeats <- function(nsample, b0, b1, size, q0, FDRcutoff=0.1){
  # predictDEFeats() generates count datasets based on input NB coefficients, 
  # and analyzes these with edgeR. The number of differentially expressed genes, 
  # the number of rejections and the number of true positives is returned
  
  if(any(is.na(c(nsample, b0, b1, size)))) stop("Input contains NA values")
  ### Check input:
  if(!all.equal(length(b0), length(b1), length(size)) | 
       (!is.null(q0) & length(b0)!=length(q0))){
    stop("Coefficient vectors are not of the same length")
  }
  require(edgeR)
  p <- length(b0)
  
  # round to nearest smaller even number:
  while(nsample%%2) nsample <- ceiling(nsample - 1)  
  if(is.null(q0)) q0 <- numeric(p)  # when likelihood is nb, not zi-nb
  
  ### Generate count matrix:  
  counts <- matrix(0, p, nsample)  # matrix of counts
  counts[, 1:(nsample/2)] <- t(sapply(1:p, function(i) 
    createZeroes(rnbinom(n=nsample/2, size=exp(size[i]), mu=exp(b0[i])), 
                 nzero=round(nsample*q0[i]/2))))
  counts[, ((nsample/2)+1):nsample] <- t(sapply(1:p, function(i) 
    createZeroes(rnbinom(n=nsample/2, size=exp(size[i]), mu=exp(b0[i]+b1[i])), 
                 nzero=round(nsample*q0[i]/2))))
  
  ### Exact testing:
  res <- DGEList(counts, norm.factors=calcNormFactors(counts, method="TMM"), 
                 group=factor(rep(1:2, each=ncol(counts)/2)))
  res <- estimateCommonDisp(res)
  res <- estimateTagwiseDisp(res)
  
  res <- exactTest(res)
  res <- topTags(res, n=nrow(counts), sort.by="none")$table
  
  
  ### Output
  return(data.frame(DE.feats=sum(b1!=0), 
                    detected.feats=sum(res$FDR <= FDRcutoff), 
                    matching.feats=sum(is.element(which(res$FDR <= FDRcutoff), 
                                                  which(b1!=0)))))
}

findTags3 <- function(nsample, b0, b1, size, q0=NULL, b1thres=0, 
                      fixedcp0=TRUE, fixedseed=TRUE, alpha=0.1){
  # findTags3() returns the number of rejections detected by exact testing
  # on a simulated dataset with input (ZI-)NB parameters, 
  # at FDR cutoff level alfa
  # nsample = number of subjects of dataset
  # b0, b1, size, q0 = ZI-NB parameters, as obtained by ShrinkBayes-analysis
  # b1thres = effect size threshold: smaller beta1s are 
  #    set to zero before data generation
  # fixedseed = boolean. fix the random seed, so that output will be 
  #             identical on consecutive calls
  if(any(is.na(c(nsample, b0, b1, size)))) stop("Input contains NA values")
  if(!all.equal(length(b0), length(b1), length(size)) | 
       (!is.null(q0) & length(b0)!=length(q0))){
    stop("Coefficient vectors are not of the same length")
  }
  require(edgeR)
  p <- length(b0)
  b1[abs(b1)<= b1thres] <- 0
  
  # round to nearest smaller even number: 
  while(nsample%%2) nsample <- ceiling(nsample - 1)
  
  fctr <- factor(rep(1:2, each=nsample/2))
  wl1 <- which(fctr=="1")  # extra overhead, but will work for factors with
  wl2 <- which(fctr=="2")  # different order
  
  if(fixedseed){
    if(!exists(".Random.seed", mode="integer", envir=globalenv())) sample(NA)
    oldSeed <- .Random.seed
    set.seed(2909)  # no visible effect on distribution
  }
  
  # Generate counts:
  sd <- matrix(0, p, nsample)
  sd[, wl1] <- t(sapply(1:p, function(i) 
    rnbinom(n=length(wl1), size=exp(size[i]), mu=exp(b0[i]))))
  sd[, wl2] <- t(sapply(1:p, function(i) 
    rnbinom(n=length(wl2), size=exp(size[i]), mu=exp(b0[i]+b1[i]))))
  
  # Create zeroes:
  if(!is.null(q0)){
    if(fixedcp0){
      for(i in 1:p) sd[i, sample(1:nsample, size=round(nsample*q0[i]))] <- 0
    }else{
      sd <- sapply(1:nsample, function(col) (runif(p) > q0)*sd[, col])
    }
  }
  
  if(fixedseed) assign(".Random.seed", oldSeed, envir=globalenv())
  
  ######## edgeR inference:
  sd <- DGEList(sd, norm.factors=calcNormFactors(sd, method="TMM"), 
                group=fctr, genes=rownames(sd))
  sd <- estimateCommonDisp(sd)
  sd <- estimateTagwiseDisp(sd)
  
  sd <- exactTest(sd)
  sd <- topTags(sd, n=nrow(sd), sort.by="none")$table
  
  return(data.frame(ActualPos=sum(b1!=0), Rejections=sum(sd$FDR <= alpha), 
                      TruePos=sum(which(b1!=0) %in% which(sd$FDR <= alpha))))
  
}


#############################
### edgeR wrapper
callEdgeR <- function(countmat, groupfact=NA, designmat=NA, normalize=FALSE, 
                      outputtable=TRUE, FDRcutoff=0.1){
  # callEdgeR() runs an edgeR analysis with tagwise dispersion estimation, 
  # from input count matrix and 
  # groupfact = two-level factor for comparison of two groups, or
  # designmat = model matrix of more complicated experimental design
  #   (using edgeR's likelihood ratio tests instead of exact tests)
  # set normalize = TRUE implements normalization by "TMM" method
  # set outputtable = FALSE to return (number of FDRs <= FDRcutoff) 
  #   instead of output table as given by topTags()-function
  if(is.na(groupfact) && is.na(designmat)[1]){
    stop("Specify experimental design \n")
  }
  if(!is.na(groupfact) && !is.na(designmat)[1]){
    cat("Multiple experimental designs specified: Higher-level design 
        from model matrix is used within GLM-producedure \n")
  }
  require(edgeR)
  
  if(normalize) countmat <- edgeRnormalization(countmat)
  
  if(!is.na(designmat[1])){  # when design is specified via model matrix
    ll <- DGEList(countmat, group=designmat[, ncol(designmat)])
    ll <- estimateGLMTrendedDisp(ll, design=designmat)
    ll <- estimateGLMTagwiseDisp(ll, design=designmat)
    
    ll <- glmFit(ll, design=designmat)
    ll <- glmLRT(ll)
    ll <- topTags(ll, n=nrow(countmat), sort.by="none")$table
  }else{  # when design is specified via factor
    ll <- DGEList(countmat, group=groupfact)
      
    ll <- estimateCommonDisp(ll) 
    ll <- estimateTagwiseDisp(ll) 
      
    ll <- exactTest(ll)
    ll <- topTags(ll, sort.by="none", n=nrow(countmat))$table
  }
  
  if(outputtable) return(ll) else return(sum(ll$FDR <= FDRcutoff))
}




