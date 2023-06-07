# Title: BRCA2 Statistical metrics
# Date: 16/05/2022

# Libraries ----
library(cvms)
library(tidyverse)
library(magrittr)
library(tibble)

# Load the algorithm output files
## The output files contain two columns: Actual and Predicted

circd_output <- c(
  "provean_circd.txt",
  "pmut_circd.txt",
  "sift_circd.txt",
  "pantherpsep_circd.txt",
  "snpsgo_circd.txt",
  "phdsnp_circd.txt",
  "predictsnp_circd.txt",
  "humdiv_circd.txt",
  "humvar_circd.txt",
  "metasnp_circd.txt"
)

data_list <- list()

for (i in seq_along(circd_output)) {
  data_list[[i]] <- read.table(circd_output[i], header = TRUE)
}

for (i in seq_along(data_list)) {
  circd_output <- basename(circd_output[i])
  variable_name <- sub(".txt$", "", circd_output)
  assign(variable_name, data_list[[i]])
}

list(circd_output)

# Generate statistical metrics ----

dataframe_names <- c(
  "provean_circd",
  "pmut_circd",
  "sift_circd",
  "pantherpsep_circd",
  "snpsgo_circd",
  "phdsnp_circd",
  "predictsnp_circd",
  "humdiv_circd",
  "humvar_circd",
  "metasnp_circd"
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
