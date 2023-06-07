# Title: BRCA2 Statistical metrics
# Date: 16/05/2022

# Libraries ----
library(cvms)
library(tidyverse)
library(magrittr)
library(tibble)

# Load the algorithm output files
## The output files contain two columns: Actual and Predicted

brca2_output <- c(
  "provean_brca2.txt",
  "pmut_brca2.txt",
  "sift_brca2.txt",
  "pantherpsep_brca2.txt",
  "snpsgo_brca2.txt",
  "phdsnp_brca2.txt",
  "predictsnp_brca2.txt",
  "humdiv_brca2.txt",
  "humvar_brca2.txt",
  "metasnp_brca2.txt"
)

data_list <- list()

for (i in seq_along(brca2_output)) {
  data_list[[i]] <- read.table(brca2_output[i], header = TRUE)
}

for (i in seq_along(data_list)) {
  brca2_output <- basename(brca2_output[i])
  variable_name <- sub(".txt$", "", brca2_output)
  assign(variable_name, data_list[[i]])
}

list(brca2_output)

# Generate statistical metrics ----

dataframe_names <- c(
  "provean_brca2",
  "pmut_brca2",
  "sift_brca2",
  "pantherpsep_brca2",
  "snpsgo_brca2",
  "phdsnp_brca2",
  "predictsnp_brca2",
  "humdiv_brca2",
  "humvar_brca2",
  "metasnp_brca2"
)

evaluated_outputs <- list()

evaluated_outputs <- lapply(dataframe_names, function(name) {
  df <- get(name)
  
  evaluated_output <- evaluate(
    df,
    target_col = "Actual",
    prediction_cols = "Predicted",
    type = "binomial"
  )
  
  evaluated_name <- paste0(name, "_evaluated")
  
  assign(evaluated_name, evaluated_output, envir = .GlobalEnv)
  
  evaluated_output
})
