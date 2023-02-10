### Simulate a data set from the true parameter

# Load source
source('scripts/utils.R')
library(Directional)

# Fix the gamma parameters
gamma_1 <- c(1, 0, 0)
gamma_2 <- c(0, 1, 0)
gamma_3 <- c(0, 0, 1)

# Parameter
k <- 2
b <- 0.5

# Random sample from Kent distribution
rkent2 <- function(n, k, m, b, gamma) {
  A <- diag(c(-b, 0, b))
  x <- Directional::rfb(n, k, m, A)
}

# Simulate
set.seed(2022)
y <- rkent2(1e2, k, gamma_1, b)
summarystats1 <- sum(y[,1] * gamma_1[1])
summarystats2 <- sum((y[,2] * gamma_2[2])^2 - (y[,3] * gamma_3[3])^2)
summary_stats <- c(summarystats1, summarystats2)
save(y, summary_stats, file = 'results/kent_simulated_data.rda')
