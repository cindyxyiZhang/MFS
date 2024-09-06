## Bonferroni
Bonferroni <- function(p_values, alpha){
  test_results <- numeric(length(p_values))
  p_vals_adjust <- p.adjust(p_values, method="bonferroni")
  test_results[which(p_vals_adjust<alpha)] <- 1
  return (test_results)
}

## Holm-Bonferroni
Holm_Bonferroni <- function(p_values, alpha){
  test_results <- numeric(length(p_values))
  p_vals_adjust <- p.adjust(p_values, method="holm")
  test_results[which(p_vals_adjust<alpha)] <- 1
  return (test_results)
}

## Oracle method
## Order of p-values known as priori 
## (mu of true alternatives not equal; so there is an order)
## Step-down: 
## p-val: smallest to largest
## Stop the 1st time that fail to reject; fail to reject this and all afterwards
## Rejection set: all previous ones


## p-val: 1,...,M (index)
## true_order: true order of p-val from smallest to largest (ordered index)
Oracle_Torder <- function(p_values, true_order, alpha){
  test_results <- rep(1, length(p_values))
  p_val_order <- p_values[true_order]
  
  kh <- which(p_val_order>=alpha) ## ordered p-val >= alpha
  if (length(kh)>0){
    kh <- min(kh)
    test_results[true_order[kh:length(true_order)]] <- 0
  }
  return (test_results)
}



## Assign label in a cyclic manner: 
## e.g. 10 elements and 3 groups
## 1->G1, 2->G2, 3->G3; 4->G1, 5->G2, 6->G3; 7->G1, 8->G2, 9->G3; 10->G1
assign_groups <- function(vec, K) {
  n <- length(vec)
  # Assign group labels in a round-robin fashion
  group_labels <- (1:n) %% K
  group_labels[group_labels == 0] <- K  # Replace 0 with K for the Kth group
  return(group_labels)
}


external_groupSD <- function(p_values, p_values_ext, ngrp, alpha){
  ## Order p-vals from external set; save indices in the original seq
  idx_ordered <- order(p_values_ext) ## smallest to largest
  ## Assign order p-vals to ngrp groups 
  group_labels <- assign_groups(idx_ordered, ngrp) ## group labels (p-val)
  
  ## Group p_vals from internal set by group labels
  ## Perform step-down procedure in each group
  test_results <- rep(1, length(p_values))
  nullSet <- lapply(1:ngrp, function(g){
    ## Extract p-vals from group g
    IndexSet_g <- idx_ordered[group_labels == g]
    pvals_g <- p_values[IndexSet_g]
    ## Hypothesis fail to reject
    
    kh_g = which(pvals_g>=alpha/ngrp)
    if (length(kh_g)>0){
      kh_g <- min(kh_g)
      IndexSet_g[kh_g:length(IndexSet_g)]
    } else {
      NULL
    }
  })
  nullSet <- unlist(nullSet)
  test_results[nullSet] <- 0
  return (test_results)
}



# external_groupSD2 <- function(p_values, p_values_ext, ngrp, alpha){
#   ## Order p-vals from external set; save indices in the original seq
#   idx_ordered <- order(p_values_ext) ## smallest to largest
#   ## Assign order p-vals to ngrp groups 
#   group_labels <- assign_groups(idx_ordered, ngrp) ## group labels (p-val)
#   
#   ## Group p_vals from internal set by group labels
#   ## Perform step-down procedure in each group
#   test_results <- rep(1, length(p_values))
#   nullSet <- lapply(1:ngrp, function(g){
#     ## Extract p-vals from group g
#     IndexSet_g <- idx_ordered[group_labels == g]
#     pvals_g <- p_values[IndexSet_g]
#     ## Hypothesis fail to reject
#     
#     kh_g = which(pvals_g>=alpha/ngrp)
#     if (length(kh_g)>0){
#       IndexSet_g[kh_g]
#     } else {
#       NULL
#     }
#   })
#   nullSet <- unlist(nullSet)
#   test_results[nullSet] <- 0
#   return (test_results)
# }








