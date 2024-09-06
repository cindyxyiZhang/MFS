rm(list = ls())
library(MASS)
library(matrixcalc) ## Efficient way generating Toeplitz matrix
library(abind)
source("Helpers_Gdk.R")

args <- commandArgs(TRUE)
n <- as.numeric(args[[1]])
run_chunk <- as.numeric(args[[2]]) 

ngrp <- 20
ntrue <- 50 # number of true alternatives

alpha = 0.05
Signals = rep(0,n) ## indicator of null and non-null
Signals[1:ntrue] <- 1 ## Put all true non-nulls at the start

cor_type <- "Toeplitz" # c("Independent", "EqualCor", "Toeplitz")
rho <- 0.9
tau <- 0.5

## mean shift vector
mu <- c(seq(3,5,length.out=ntrue), numeric(n-ntrue))

niter = 10
FWER = array(0,c(4,1,niter))
power = array(0,c(4,1,niter))


for(iter in 1:niter){
  set.seed(2024+niter*(run_chunk-1)+iter)
  
  if (cor_type=="Independent"){
    P = 1-pnorm(rnorm(n) + mu*Signals)
    Pext = 1-pnorm(rnorm(n) + mu*Signals) ## p-values from external set
  } else if (cor_type=="EqualCor"){
    Sigma <- matrix(rho, n, n)
    diag(Sigma) <- 1  
    P = 1-pnorm(as.vector(rnorm(n)%*% chol(Sigma)) + mu*Signals)
    Pext = 1-pnorm(as.vector(rnorm(n)%*% chol(Sigma)) + mu*Signals)
  } else if (cor_type=="Toeplitz"){
    first_row <- tau^(0:(n-1))
    # Generate the Toeplitz covariance matrix
    Sigma <- toeplitz(first_row)
    P = 1-pnorm(as.vector(rnorm(n)%*% chol(Sigma)) + mu*Signals)
    Pext = 1-pnorm(as.vector(rnorm(n)%*% chol(Sigma)) + mu*Signals)
  }
  
  result_Bonferroni <- Bonferroni(P, alpha)
  result_Holm <- Holm_Bonferroni(P, alpha)
  result_Oracle <- Oracle_Torder(P, c(50:1, 51:n), alpha)
  result_extSD <- external_groupSD(P, Pext, ngrp, alpha)
  
  results=cbind(result_Bonferroni, result_Holm, result_Oracle, result_extSD)
  
  for(method in 1:4){	
    FWER[method,1,iter]=ifelse(sum(results[,method]*(1-Signals)) > 0, 1, 0)
    
    power[method,1,iter]=sum(results[,method]*Signals)/sum(Signals)
    
  }
}

save(FWER, power, file=paste(c("Output/simulation_Gdk_Toeplitz_p", n, 
                               "_", run_chunk, ".rda"), collapse = ""))

