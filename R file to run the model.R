########################################################################
#### R code file to run the spatial misalignment model on simulated data
########################################################################

set.seed(1)

#####################################
#### Load the libraries and functions
#####################################
#### Libraries
library(sp)
library(rgdal)
library(truncdist)
library(coda)
library(truncnorm)
library(Rcpp)
library(MASS)
library(MCMCpack)


#### Functions
source('print.CARBayes.R')
Rcpp::sourceCpp('CARBayes.cpp')
source('common.functions.R')
source('poisson.MVlerouxCARdasmoothalpha.R')



###############################
#### Read in the simulated data
###############################
load(file="Data GPlevel.Rdata")
load(file="Data IZlevel.Rdata")
load(file="GPIZ populations.Rdata")
K.GP <- nrow(gp.dat)
K.IZ <- nrow(IZ.dat)



##################################
#### Read in the IZ level W matrix
##################################
load(file="WIZ.Rdata")



################################################################
#### Compute the population weights and IZ level expected counts
################################################################
weights.e <- GPIZ / matrix(rep(rowSums(GPIZ), K.IZ), nrow=K.GP, ncol=K.IZ, byrow=FALSE)
e.IZ.gp <- t(weights.e) %*%  gp.dat$Egp



#####################################################
#### Compute the initial Y.IZ  counts for the GP data
#####################################################
Y.intersections <- array(0, c(K.GP, K.IZ))
  for(s in 1:K.GP)
  {
  probs.unscaled <- weights.e[s, ]
  probs <- probs.unscaled / sum(probs.unscaled)
  Y.intersections[s, ] <- as.numeric(rmultinom(n=1, size=gp.dat$Ygp[s], prob=probs))
  }
Y.IZ.gp <- round(apply(Y.intersections,2,sum))   



####################
#### MCMC quantities
####################
burnin <- 100000
n.sample <- 1100000
thin <- 200
n.da <- 20



#######################
#### Fit the model data
#######################
#### Specify the variables and formula
Y.final <- cbind(Y.IZ.gp, IZ.dat$Yhosp, IZ.dat$Ydeath)
e.final <- cbind(e.IZ.gp, IZ.dat$Ehosp, IZ.dat$Edeath)
simd <- IZ.dat$EmpRate
form <- Y.final ~ offset(log(e.final)) + simd
Y.region <- cbind(gp.dat$Ygp, rep(NA, K.GP), rep(NA, K.GP))


#### Run the model
chain1 <- poisson.MVlerouxCARdasmoothalpha(formula=form, data=NULL, W=W, Y.region=Y.region, weights=weights.e, alpha=0.9,  burnin=burnin, n.sample=n.sample, thin=thin, n.da=n.da)
chain2 <- poisson.MVlerouxCARdasmoothalpha(formula=form, data=NULL, W=W, Y.region=Y.region, weights=weights.e, alpha=0.9,  burnin=burnin, n.sample=n.sample, thin=thin, n.da=n.da)
chain3 <- poisson.MVlerouxCARdasmoothalpha(formula=form, data=NULL, W=W, Y.region=Y.region, weights=weights.e, alpha=0.9,  burnin=burnin, n.sample=n.sample, thin=thin, n.da=n.da)


#### Print the results
print(chain1)
print(chain2)
print(chain3)


#### Model fit
chain1$modelfit
chain2$modelfit
chain3$modelfit