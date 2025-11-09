# library(MCMCpack)
# library(MCMCglmm)
# #coda
# #parallel

# data(BTdata)

# devtools::document("/Users/joelpick/github/mcMCMCglmm")

# prior=list(
# 	R=list(
# 		R1=list(V = 1, nu = 0.002)
# 	),
# 	G=list(
# 		G1=list(V = 1, nu = 0.002),
# 		G2=list(V = 1, nu = 0.002)
# 	)
# )
# mod <- mcMCMCglmm(tarsus~back, random=~dam+fosternest, data=BTdata,prior=prior)

# # prior2=list(
# # 	R=list(
# # 		R1=list(V = diag(2), nu = 2.002)
# # 	),
# # 	G=list(
# # 		G1=list(V = diag(2), nu = 2.002),
# # 		G2=list(V = diag(2), nu = 2.002)
# # 	)
# # )
# # mod <- MCMCglmm(cbind(tarsus,back)~1, random=~us(trait):dam+us(trait):fosternest,rcov=~us(trait):units,data=BTdata, family=c("gaussian","gaussian"),prior=prior2)


# ##how to deal with covariances 
# ## something to indicate whether they are a variance or covariance, and then plot the priors for the variances
# summary(mod)
# plot(mod)

# # post<-mod$VCV[[1]]
# # x=seq(0.005,max(post),length.out=100)
# # prior_density<-do.call(dprior,c(list(x=x),prior$G$G1))

# # pp_plot(x,prior_density,post)
# ## try  rIW() from MCMCglmm

# which_diag <- function(n){
# 	x<-diag(n)
# 	as.vector(!(lower.tri(x)|upper.tri(x)))
# }

# var_cols<- c(sapply(c(mod$Random$nfl,mod$Residual$nfl),which_diag), recursive=TRUE)

# var_post <- do.call(rbind,mod$VCV[,var_cols,])

# x <- seq(min(var_post),max(var_post),length.out=100)
# g_priors<-do.call(cbind,lapply(mod$prior$G,function(prior){
# 	sapply(1:sqrt(length(prior$V)), function(variable){
# 		do.call(dprior,c(list(x=x,variable=variable),prior))
# 	})
# }))
# r_priors<-do.call(cbind,lapply(mod$prior$R,function(prior){
# 	sapply(1:sqrt(length(prior$V)), function(variable){
# 		do.call(dprior,c(list(x=x,variable=variable),prior))
# 	})
# }))
# var_prior <- cbind(g_priors,r_priors)

# ncol(var_prior)
# ncol(var_post)

# par(mfrow=c(2,3))
# for(i in 1:ncol(var_prior)) pp_plot(x,prior=var_prior[,i],posterior=var_post[,i])

# prior3=list(
# 	R=list(
# 		R1=list(V = diag(2), nu = 1.002)
# 	),
# 	G=list(
# 		G1=list(V = diag(c(1,2)), nu = 1.002),
# 		G2=list(V = diag(c(1,2)), nu = 1.1)
# 	)
# )


# post<-mod$VCV[[1]]
# x=seq(0.005,max(post),length.out=100)
# prior_density<-do.call(dprior,c(list(x=x),prior$G$G1))


