# Title: BRCA2 Statistical metrics
# Date: 16/05/2022

# Libraries ----
library(cvms)
library(tidyverse)
library(magrittr)
library(tibble)

# Load the algorithm output files
## The output files contain two columns: Actual and Predicted

ep_output <- c(
  "provean_ep.txt",
  "pmut_ep.txt",
  "sift_ep.txt",
  "pantherpsep_ep.txt",
  "snpsgo_ep.txt",
  "phdsnp_ep.txt",
  "predictsnp_ep.txt",
  "humdiv_ep.txt",
  "humvar_ep.txt",
  "metasnp_ep.txt"
)

data_list <- list()

for (i in seq_along(ep_output)) {
  data_list[[i]] <- read.table(ep_output[i], header = TRUE)
}

for (i in seq_along(data_list)) {
  ep_output <- basename(ep_output[i])
  variable_name <- sub(".txt$", "", ep_output)
  assign(variable_name, data_list[[i]])
}

list(ep_output)

# Generate statistical metrics ----

dataframe_names <- c(
  "provean_ep",
  "pmut_ep",
  "sift_ep",
  "pantherpsep_ep",
  "snpsgo_ep",
  "phdsnp_ep",
  "predictsnp_ep",
  "humdiv_ep",
  "humvar_ep",
  "metasnp_ep"
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
