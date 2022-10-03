### Simulate a data set from the true parameter

# Load source
source('scripts/utils.R')

# Initialize parameters
n <- 1000
d <- 2
k <- 5

# Setup
dims <- rep(n,d)
if (d == 2) {
  mask <- matrix(1,dims[1],dims[2])
  neigh <- getNeighbors(mask, c(2,2,0,0))
  edges <- getEdges(mask, c(2,2,0,0))
  bcrit <- log(1 + sqrt(k)) # Potts (1952)
}
block <- getBlocks(mask, 2)
maxS <- nrow(edges)

# Simulate an image
set.seed(8)
beta_true <- rnorm(1, bcrit, 0.05)
for (i in 1:5){
set.seed(i)
sim <- swNoData(beta_true, k, neigh, block, 501)
summary_stats <- sim$sum[501]
save(summary_stats, file = paste0('results/potts_simulated_data_', i, '.rda'))
}