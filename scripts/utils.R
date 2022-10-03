### Load the required packages and functions

# reticulate::use_condaenv("TFv1-15", required=TRUE)

library(bayesImageS)
library(coda)
library(devtools)
library(dplyr)
library(doParallel)
library(GiRaF)
library(ggplot2)
library(gridExtra)
# library(PAWL)
load_all('deepsurrogate')

## Function to count neighboring pixels
CountNeighbors <- function(matrix){
  neigh <- 0
  for (i in 1:nrow(matrix)){
    for (j in 1:ncol(matrix)){
      if (j > 1){
        if (matrix[i,j] == matrix[i, j-1]){neigh <- neigh + 1}
      }
      if (i > 1){
        if (matrix[i,j] == matrix[i-1, j]){neigh <- neigh + 1}
      }
    }
  }
  return(neigh)
}
