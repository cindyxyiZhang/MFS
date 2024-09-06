#!/bin/bash

# Array of s*level2 folders
s_levels=("s2level2") # "s4level2"

# Array of cope*.feat folders
cope_levels=("cope4.feat")
# "cope9.feat" "cope10.feat" "cope11.feat" "cope12.feat" "cope13.feat"

# Maximum number of parallel jobs
MAX_JOBS=2

# Loop through each combination and call the R script
for s in "${s_levels[@]}"; do
  for cope in "${cope_levels[@]}"; do
    ((i=i%MAX_JOBS)); ((i++==0)) && wait
    Rscript groupLevel_pval_fewSubj_BonfAdd.R $s $cope &
  done
done

wait

