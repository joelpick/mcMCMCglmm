#mcMCMCglmm

#' @title mcMCMCglmm
#' @description Function to run multi-chain MCMCglmm models in parallel
#' @param ... arguments passed to MCMCglmm
#' @param chains number of chains
#' @param n.cores number of cores
#' @param prior prior argument to MCMCglmm
#' @param verbose verbose argument to MCMCglmm
#' @details Pass the same arguements to MCMCglmm. It a smilar MCMCglmm is returned, the only difference is the slots taht previously had single posteriors for each parameter, have mcmc.list objects with multiple posteriors for each parameter (e.g. Sol and VCV)
#' @author Joel Pick - joel.l.pick@gmail.com
#' @export
mcMCMCglmm <- function(..., chains=4, n.cores=1, prior = NULL, verbose=FALSE){
	mods <- parallel::mclapply(1:chains,function(i){
		mod <- MCMCglmm::MCMCglmm(..., prior = prior, verbose=verbose)
		cat("-- Chain", i, "complete --\n")
		return(mod)
	}, mc.cores=n.cores)

	## combine models together to get mcmc.lists for each component with a 
	mod_return <- combine_models(mods)

	## save the prior, and then can plot it.
	mod_return$prior <- prior	
	# from the help file:
	# The defaults are ‘nu’=0, ‘V’=1,‘alpha.mu’=0, and ‘alpha.V’=0. When ‘alpha.V’ is non-zero, parameter expanded algorithms are used.

	class(mod_return) <- c("mcMCMCglmm","MCMCglmm")
	return(mod_return)
}


## internal function that combines a list with MCMCglmm model outputs run with same call to MCMCglmm, to give one object, with components that are posteriors grouped into mcmc.list objects
combine_models <- function(model_list){

	mod_combine <- model_list[[1]]
	mod_combine$Sol<- coda::mcmc.list(lapply(model_list,function(model_list) model_list$Sol))
	if(!is.null(model_list[[1]]$Lambda)){
		mod_combine$Lambda<- coda::mcmc.list(lapply(model_list,function(model_list) model_list$Lambda))
	}
	mod_combine$VCV<- coda::mcmc.list(lapply(model_list,function(model_list) model_list$VCV))
	if(!is.null(model_list[[1]]$CP)){
		mod_combine$CP<- coda::mcmc.list(lapply(model_list,function(model_list) model_list$CP))
	}
	if(!is.null(model_list[[1]]$Liab)){
		mod_combine$Liab<- coda::mcmc.list(lapply(model_list,function(model_list) model_list$Liab))
	}
	if(!is.null(model_list[[1]]$ThetaS)){
		mod_combine$ThetaS<- coda::mcmc.list(lapply(model_list,function(model_list) model_list$ThetaS))
	}
	mod_combine$Deviance <- coda::mcmc.list(lapply(model_list,function(model_list) model_list$Deviance))
	mod_combine$DIC <- sapply(model_list,function(model_list) model_list$DIC)

	mod_combine
}

## Rhat from BDA2, Gelman et al., 2003
## takes matrix of one param with all chains
Rhat <- function(x){
	n_samples <- nrow(x)
 	chain_mean <- colMeans(x)
  chain_var <- apply(x,2,stats::var)
  
  var_between <- n_samples * stats::var(chain_mean)
  var_within <- mean(chain_var)
  sqrt((var_between / var_within + n_samples - 1) / n_samples)
}

## Rhat on mcmc.list
Rhat_multi <- function(post){
	post_array <- as.array(post)
	if(length(dim(post_array))==3){
		apply(post_array,2,Rhat)
	}else{
		Rhat(post_array)
	}
}

# 'p-value' from MCMC chains
pMCMC<- function(x){
	n_samples <- dim(x)[1]
	2 * pmax(
		0.5/n_samples, 
		pmin(colSums(x > 0)/n_samples, 
			1 - colSums(x > 0)/n_samples
		)
	)
}

## returns summary from mcmc.list 
results_out<- function(x, pMCMC=FALSE){
	x_single <- coda::mcmc(do.call(rbind,x))
	out <- cbind(
  	colMeans(x_single), 
    coda::HPDinterval(x_single), 
    coda::effectiveSize(x), 
    if(pMCMC){pMCMC(x_single)},
    Rhat_multi(x)
  )
  colnames(out) <- c("post.mean", "l-95% CI", "u-95% CI", "eff.samp", if(pMCMC){"pMCMC"},"Rhat")
  out
}

#' @title summary.mcMCMCglmm
#' @description Function to provide summary for mcMCMCglmm model, providing give summary stats across multiple chains. Gives the same output as summary.MCMCglmm, with the additional of Rhat
#' @param object mcMCMCglmm object
#' @param random logical: should the random effects be summarised
#' @param ... Additional arguments
#' @details Provides the same output as MCMCglmm summary function with the addition of Rhat
#' @export

summary.mcMCMCglmm <- function (object, random = FALSE, ...){

  DIC <- object$DIC
  fixed.formula <- object$Fixed$formula
  nF <- object$Fixed$nfl
  nL <- object$Fixed$nll
  if (random) {
      nF <- sum(rep(object$Random$nrl, object$Random$nfl)) + 
          nF
      if (nF != dim(object$Sol)[2]) {
          stop("random effects not saved and cannot be summarised")
      }
  }
  solutions <- results_out(object$Sol[, 1:nF, drop = FALSE], pMCMC=TRUE)
  if (nL > 0) {
    solutions <- rbind(solutions, 
			results_out(object$Lambda, pMCMC=TRUE)
    )
  }
  
  random.formula = object$Random$formula
  residual.formula = object$Residual$formula
  gterms <- sum(object$Random$nfl^2)
  rterms <- sum(object$Residual$nfl^2)
  covariances <- results_out(object$VCV)

  if (gterms > 0) {
      Gcovariances <- covariances[1:gterms, , drop = FALSE]
  } else {
      Gcovariances <- NULL
  }
  Rcovariances <- covariances[gterms + 1:rterms, , drop = FALSE]
  cstats <- attr(object$VCV[[1]], "mcpar")
  cstats[4] <- dim(object$VCV[[1]])[1]
  if (is.null(object$CP)) {
      cutpoints <- NULL
  } else {
      cutpoints <- results_out(object$CP)
  }
  if (is.null(object$ThetaS)) {
      theta_scale <- NULL
  } else {
      theta_scale <- results_out(object$ThetaS)
  }
  if (is.null(object$Random$nrt)) {
      Gterms <- NULL
  } else {
      Gterms <- rep(rep(1:length(object$Random$nrt), object$Random$nrt), object$Random$nfl^2)
  }
  
  Rterms <- rep(rep(1:length(object$Residual$nrt), object$Residual$nrt), object$Residual$nfl^2)

  output <- list(
  	DIC = DIC, 
  	fixed.formula = fixed.formula, 
    random.formula = random.formula, 
    residual.formula = residual.formula, 
    solutions = solutions, 
    Gcovariances = Gcovariances, 
    Gterms = Gterms, 
    Rcovariances = Rcovariances, 
    Rterms = Rterms, 
    cstats = cstats, 
    cutpoints = cutpoints, 
    theta_scale = theta_scale)
  attr(output, "class") <- c("summary.MCMCglmm", "list")
  output
}

