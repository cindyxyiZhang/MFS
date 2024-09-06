## Exploration: gain and loss by splitting data to estimate ordering vs. Bonferroni on full

rm(list = ls())
library(MASS)
library(matrixcalc) ## Efficient way generating Toeplitz matrix
library(abind)
source("Helpers_Gdk.R")

args <- commandArgs(TRUE)
ngrp <- as.numeric(args[[1]]) # group no. hyperparameter: 10, 50, 100, 500, 1000 
Nsubj <- as.numeric(args[[2]]) # sample size 
n <- as.numeric(args[[3]]) ## number of hypothesis
ntrue <- 50 # number of true alternatives
split_prop <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
# 0.6 0.7 0.8

alpha = 0.05
Signals = rep(0,n) ## indicator of null and non-null
Signals[1:ntrue] <- 1 ## Put all true non-nulls at the start

cor_type <- "Independent" # c("Independent", "EqualCor", "Toeplitz")
rho <- 0.9
tau <- 0.5

## mean shift vector
mu <- c(seq(0.2,1,length.out=ntrue), numeric(n-ntrue))

## Save output for each ngrp: 4 methods in comparison 
## Bind results by 2nd dimension across ngrp values
niter = 1000
FWER = array(0, c((length(split_prop)+2), 1, niter))
power = array(0, c((length(split_prop)+2), 1, niter))
dat_sd <- 1
for(iter in 1:niter){
  set.seed(2024+iter)
  if (cor_type=="Independent"){
    ## Create data matrix
    X <-  do.call("rbind", replicate(Nsubj, (rnorm(n, sd=dat_sd) + mu*Signals), simplify = FALSE))
    # X <-  do.call("rbind", replicate(Nsubj, (rt(n, df=3) + mu*Signals), simplify = FALSE))
  } else if (cor_type=="EqualCor"){
    Sigma <- matrix(rho, n, n)
    diag(Sigma) <- 1  
    X <-  do.call("rbind", replicate(Nsubj, (as.vector(rnorm(n, sd=dat_sd)%*% chol(Sigma)) + mu*Signals), simplify = FALSE)) 
  } else if (cor_type=="Toeplitz"){
    first_row <- tau^(0:(n-1))
    # Generate the Toeplitz covariance matrix
    Sigma <- toeplitz(first_row)
    X <-  do.call("rbind", replicate(Nsubj, (as.vector(rnorm(n, sd=dat_sd)%*% chol(Sigma)) + mu*Signals), simplify = FALSE)) 
  }
  
  ## p-val on complete data
  Pval_complete <- lapply(1:n, function(j){
    # Tj <- mean(X[,j])/(sd(X[,j])/sqrt(Nsubj))
    # p_val_j <- 2*(1-pnorm(abs(Tj)))
    
    # t test
    Tj <- t.test(X[,j], mu = 0)
    p_val_j <- Tj$p.value
    p_val_j
  })
  Pval_complete <- unlist(Pval_complete)
  
  result_Bonferroni <- Bonferroni(Pval_complete, alpha)
  result_Oracle <- Oracle_Torder(Pval_complete, c(ntrue:1, (ntrue+1):n), alpha)
  
  ## Multi-step-down via sample splitting
  result_ssMSD <- lapply(split_prop, function(prop){
    S1_index <- sample(1:Nsubj, floor(Nsubj*prop), replace = FALSE)# subjects for testing; 
    # remaining for ordering estimation
    n1 <- length(S1_index) # For testing
    P <- lapply(1:n, function(j){
      # Tj <- mean(X[S1_index, j])/(sd(X[S1_index, j])/sqrt(n1))
      # p_val_j <- 2*(1-pnorm(abs(Tj)))
      
      Tj <- t.test(X[S1_index, j], mu = 0)
      p_val_j <- Tj$p.value
      p_val_j
    })
    P <- unlist(P)
    
    Pext <- lapply(1:n, function(j){
      # Tj <- mean(X[-S1_index, j])/(sd(X[-S1_index, j])/sqrt(Nsubj-n1))
      # p_val_j <- 2*(1-pnorm(abs(Tj)))
      
      Tj <- t.test(X[-S1_index, j], mu = 0)
      p_val_j <- Tj$p.value
      p_val_j
    })
    Pext <- unlist(Pext)
    external_groupSD(P, Pext, ngrp, alpha)
  })
  result_ssMSD <- do.call("cbind", result_ssMSD)
  
  results=cbind(result_Bonferroni, result_Oracle, result_ssMSD)
  for(method in 1:ncol(results)){	
    FWER[method,1,iter]=ifelse(sum(results[,method]*(1-Signals)) > 0, 1, 0)
    
    power[method,1,iter]=sum(results[,method]*Signals)/sum(Signals)
    
  }
}

save(FWER, power, file=paste(c("Output/simulationSS_n", Nsubj, "_p", n, "_G", ngrp, 
                               "_indep.rda"), collapse = ""))


# Combine results
Nsubj <- 500 #100 300 500
n <- 50000  #1000 10000 50000
powerBYngrp <- list()
fwerBYngrp <- lapply(c(1, 5, 10, 20, 30, 50, 60, 80, 100, 200, 500, 1000), function(ngrp){
  load(file=paste(c("/Volumes/XYZ_HD/MultipleTesting/Simulations/MSD_SS_Independent/simulationSS_n", Nsubj, "_p", n, "_G", ngrp,
                    "_indep.rda"), collapse = ""))
  powerBYngrp <<- c(powerBYngrp, list(power))
  FWER
})
FWER <- abind(fwerBYngrp, along = 2)
power <- abind(powerBYngrp, along = 2)
save(power, FWER, file=paste(c("/Volumes/XYZ_HD/MultipleTesting/Combine/simulationSS_n", Nsubj, "_p", n, ".rda"),
                             collapse = ""))



## Toepliz
Nsubj <- 100 #100 300 500
powerBYngrp <- list()
fwerBYngrp <- lapply(c(1, 20, 50, 60, 80, 100, 500, 1000), function(ngrp){
  load(file=paste(c("/Volumes/XYZ_HD/MultipleTesting/Simulations/MSD_SS_Toeplitz/simulationSS_n", Nsubj, "_G", ngrp,
                    "_Toeplitz.rda"), collapse = ""))
  powerBYngrp <<- c(powerBYngrp, list(power))
  FWER
})
FWER <- abind(fwerBYngrp, along = 2)
power <- abind(powerBYngrp, along = 2)
save(power, FWER, file=paste(c("/Volumes/XYZ_HD/MultipleTesting/Combine/simulationSS_n", Nsubj, "_p", n, ".rda"),
                             collapse = ""))

