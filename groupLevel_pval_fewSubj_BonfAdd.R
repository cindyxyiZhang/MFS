rm(list = ls())

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
s_level <- args[1]
cope_level <- args[2]


# Extract the numeric part from the cope_level string for file names
cope_num <- sub("cope([0-9]+)\\.feat", "\\1", cope_level)

path <- sprintf("/users/xzhang3/HCP/TaskLanguage/%s/%s", s_level, cope_level)
file_cope <- sprintf("%s/cope%sCombine.txt", path, cope_num)

# Read the entire datasets and convert to matrix
cope_data <- as.matrix(fread(file_cope))

# Check and remove columns with NA in either file, only if NAs are present
na_columns <- which(colSums(is.na(cope_data)) > 0)
if (length(na_columns) > 0) {
  cope_data <- cope_data[, -na_columns]
}

# Set seed for reproducibility for internal analysis
set.seed(2024)
internal_cols <- sample(ncol(cope_data), 20)
internal_data <- cope_data[, internal_cols]

# Save internal results
output_dir <- sprintf("%s/pvals", path)
#dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Remaining columns for external analysis
external_cols <- setdiff(1:ncol(cope_data), internal_cols)

# Define column subsets
set.seed(2023)  # New seed for external analysis
col_subsets <- list(
  "20" = sample(external_cols, 20),
  "50" = sample(external_cols, 50),
  "100" = sample(external_cols, 100),
  "300" = sample(external_cols, 300)
  # "500" = sample(external_cols, 500),
  # "869" = external_cols
)

# Perform tests on each subset row-wise and save results
for (subset_name in names(col_subsets)) {
  subset_data <- cope_data[, col_subsets[[subset_name]]]
  ## Combine with original internal data; boost size for Bonferroni
  subset_data <- cbind(internal_data, subset_data)
  
  pvals_internal_add <- apply(subset_data, 1, function(row) {
    # Ttest <- t.test(row, mu = 0, alternative = "two.sided")
    # Ttest$p.value
    wilcoxTest <- wilcox.test(row, mu = 0, alternative = "two.sided")
    wilcoxTest$p.value
  })
  output_file_internal_add <- sprintf("%s/pvalsWilcox_internalAdd_%s.rda", output_dir, subset_name)
  save(pvals_internal_add, file=output_file_internal_add)
}





