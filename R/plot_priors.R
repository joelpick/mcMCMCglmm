#' @title dprior
#' @description Density function for MCMCglmm priors 
#' @param x vector of values for which the density is returned
#' @param V V in MCMCglmm prior specification 
#' @param nu nu in MCMCglmm prior specification 
#' @param alpha.mu alpha.mu in MCMCglmm prior specification
#' @param alpha.V alpha.V in MCMCglmm prior specification
#' @param variable when multivariate priors specified, which variable they apply to 
#' @details When V and nu only are specified, gives an inverse wishart prior, and when a covariance matrix is specified for V, it gives marginal prior on variance of specified variable (p12 and p72 of course notes).
#' When alpha.mu and alpha.V are specified, it returns a parameter expanded prior (p129 of course notes).
#' 
dprior <- function(x,V,nu,alpha.mu=NULL,alpha.V=NULL,variable=1){
	if((is.null(alpha.mu)&!is.null(alpha.V))|(!is.null(alpha.mu)&is.null(alpha.V))) warning("Both alpha.mu and alpha.V have to be specified for an parameter expanded prior - inverse wishart prior returned.")

	if(is.null(alpha.mu)|is.null(alpha.V)){
		# inverse wishart prior, with covariance matrix gives marginal prior on variance of specified variable (p12 and p72 of course notes)
		if(!is.matrix(V)) V <- as.matrix(V)
		nu.ast <- nu - dim(V)[1] + 1
		if(nu.ast<0) stop("nu must be more than number of variables")
		V.ast <- V[variable,variable] * (nu/nu.ast) 
		# print(c(V.ast=V.ast,nu.ast=nu.ast))
		MCMCpack::dinvgamma(x, shape = nu.ast/2, scale = (nu.ast * V.ast)/2)
	}else{
		if(!is.matrix(alpha.V)) alpha.V <- as.matrix(alpha.V)

		# parameter expanded prior - p129
		stats::df(x/alpha.V[variable,variable], df1 = 1, df2 = nu, ncp = (alpha.mu[variable]^2)/alpha.V[variable,variable])
	}
# but see footnote p105

}


#' @title pp_plot
#' @description Function to plot prior and posterior 
#' @param x vector of values for which the density is returned
#' @param prior vector of the prior density at each of x
#' @param posterior posterior distribution with which to compare the prior
#' @param main plot title
#' @param xlab label for x axis
#' @param ylim limits for y axis
#' @param xlim limits for x axis
#' @details Where x is x sequence for the x values, and also is used as the breaks in the histogram, prior is the prior density at x (i.e. applying the prior density function to x), and posterior is your posterior distribution

pp_plot<-function(x,prior,posterior,main="",xlab="",ylim=NULL,xlim=NULL){
	scale_d <- (length(posterior))/(sum(prior))
	prior_y <- prior*scale_d
	hist_out <- graphics::hist(posterior,breaks=x, plot=FALSE)
	if(is.null(ylim)) ylim <- c(0,max(prior_y,hist_out$counts))
	if(is.null(xlim)) xlim <- range(x)
	hist_out <- graphics::hist(posterior, xlim=xlim,breaks=x, col="grey",main=main,xlab="", ylim=ylim)
	graphics::lines(prior_y~x, col="red", lwd=2)
	graphics::mtext(xlab, 1, line=3, cex=1.5)
}


# Example:	
	
# pp_plot(
# 	x = seq(4,40,length.out=101),
# 	prior = dunif(seq(4,40,length.out=101), 4,40),
# 	posterior=model_zpost[,"nu_nest”]
# )