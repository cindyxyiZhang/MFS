rm(list = ls())
library(MASS)
library(matrixcalc) ## Efficient way generating Toeplitz matrix
library(abind)
source("Helpers_Gdk.R")


args <- commandArgs(TRUE)
sd_ext <- as.numeric(args[[1]]) 
n <- 1000 ## number of hypothesis
ntrue <- 50 # number of true alternatives

alpha = 0.05
Signals = rep(0,n) ## indicator of null and non-null
Signals[1:ntrue] <- 1 ## Put all true non-nulls at the start

cor_type <- "Toeplitz" # c("Independent", "EqualCor", "Toeplitz")
rho <- 0.9
tau <- 0.5

## mean shift vector
mu <- c(seq(3,5,length.out=ntrue), numeric(n-ntrue))

ngrpSeq <- c(1, 5, 10, 20, 30, 50, 60, 80, 100, 500, 1000)

## Save output for each ngrp: 2nd dimension--methods
## Bind results by 2nd dimension across ngrp values
niter = 1000
FWER = array(0, c(length(ngrpSeq),1,niter))
power = array(0, c(length(ngrpSeq),1,niter))

for(iter in 1:niter){
  set.seed(2024+iter)
  if (cor_type=="Independent"){
    P = 1-pnorm(rnorm(n) + mu*Signals)
    Pext = 1-pnorm(rnorm(n, mean=0, sd=sd_ext) + mu*Signals) ## p-values from external set
  } else if (cor_type=="EqualCor"){
    Sigma <- matrix(rho, n, n)
    diag(Sigma) <- 1  
    P = 1-pnorm(as.vector(rnorm(n)%*% chol(Sigma)) + mu*Signals)
    Pext = 1-pnorm(as.vector(rnorm(n, mean=0, sd=sd_ext)%*% chol(Sigma)) + mu*Signals)
  } else if (cor_type=="Toeplitz"){
    first_row <- tau^(0:(n-1))
    # Generate the Toeplitz covariance matrix
    Sigma <- toeplitz(first_row)
    P = 1-pnorm(as.vector(rnorm(n)%*% chol(Sigma)) + mu*Signals)
    Pext = 1-pnorm(as.vector(rnorm(n, mean=0, sd=sd_ext)%*% chol(Sigma)) + mu*Signals)
  }
  
  results <- lapply(ngrpSeq, function(ngrp){
    external_groupSD(P, Pext, ngrp, alpha)
  })
  results <- do.call("cbind", results)

  for(method in 1:length(ngrpSeq)){	
    FWER[method,1,iter]=ifelse(sum(results[,method]*(1-Signals)) > 0, 1, 0)
    
    power[method,1,iter]=sum(results[,method]*Signals)/sum(Signals)
    
  }
}

# MBSE: Mixture(M) of Bonferroni(B) and Step-down(S) procedure with External(E) data
save(FWER, power, file=paste(c("Output/simulation_GdkMBSE_", "sd_", sd_ext,
                               "_Toeplitz.rda"), collapse = ""))


# ## Combine results
powerBYngrp <- list()
fwerBYngrp <- lapply(c(0.2, 0.5, 0.8, 1, 1.5, 2, 3, 4), function(sd_ext){
  load(file=paste(c("/Volumes/XYZ_HD/MultipleTesting/simulation_Gdk_", 50, "sd_", sd_ext,
                    "_Toeplitz.rda"), collapse = ""))
  powerBYngrp <<- c(powerBYngrp, list(power))
  FWER
})
FWER <- abind(fwerBYngrp, along = 2)
power <- abind(powerBYngrp, along = 2)
## Only extract results from Bonferroni, Holm-Bonferroni and Oracle
## Does NOT depend on group and external data
## Extract the 1st three rows 
FWER  <- FWER[1:3, , ]
power <- power[1:3, , ]

## Results from various group sizes and external error variance
powerBYngrp <- list()
fwerBYngrp <- lapply(c(0.2, 0.5, 0.8, 1, 1.5, 2, 3, 4), function(sd_ext){
  load(file=paste(c("/Volumes/XYZ_HD/MultipleTesting/simulation_GdkMBSE_", "sd_", sd_ext,
                    "_Toeplitz.rda"), collapse = ""))
  powerBYngrp <<- c(powerBYngrp, list(power))
  FWER
})
FWER_gpExtError <- abind(fwerBYngrp, along = 2)
power_gpExtError <- abind(powerBYngrp, along = 2)

## Combine with other methods
FWER <- abind(FWER, FWER_gpExtError, along = 1) # 1st dim=14; 2nd dim=8; 3rd dim=1000
power <- abind(power, power_gpExtError, along = 1)

save(power, FWER, file="/Volumes/XYZ_HD/MultipleTesting/Combine/simulation_Gdk_ToeplitzErrorGroupVary.rda")





