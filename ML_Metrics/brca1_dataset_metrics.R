# Title: BRCA2 Statistical metrics
# Date: 16/05/2022

# Libraries ----
library(cvms)
library(tidyverse)
library(magrittr)
library(tibble)

# Load the algorithm output files
## The output files contain two columns: Actual and Predicted

brca1_output <- c(
  "provean_brca1.txt",
  "pmut_brca1.txt",
  "sift_brca1.txt",
  "pantherpsep_brca1.txt",
  "snpsgo_brca1.txt",
  "phdsnp_brca1.txt",
  "predictsnp_brca1.txt",
  "humdiv_brca1.txt",
  "humvar_brca1.txt",
  "metasnp_brca1.txt"
)

data_list <- list()

for (i in seq_along(brca1_output)) {
  data_list[[i]] <- read.table(brca1_output[i], header = TRUE)
}

for (i in seq_along(data_list)) {
  brca1_output <- basename(brca1_output[i])
  variable_name <- sub(".txt$", "", brca1_output)
  assign(variable_name, data_list[[i]])
}

list(brca1_output)

# Generate statistical metrics ----

dataframe_names <- c(
  "provean_brca1",
  "pmut_brca1",
  "sift_brca1",
  "pantherpsep_brca1",
  "snpsgo_brca1",
  "phdsnp_brca1",
  "predictsnp_brca1",
  "humdiv_brca1",
  "humvar_brca1",
  "metasnp_brca1"
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
