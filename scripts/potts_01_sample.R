### Simulate the pseudo-data to fit the surrogate model

# Load source
source('scripts/utils.R')
load('results/potts_simulated_data.rda')
registerDoParallel(cores = 64L)

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

############ Sample from a grid

t1 <- proc.time()
beta <- seq(0.9, 1.3, by = 0.004)
matu <- matrix(0, 1000, length(beta))
sample_data <- foreach (i = 1:length(beta)) %dopar% {
  set.seed(i)
  res2 <- swNoData(beta[i], k, neigh, block, 1000)
  matu[,i] <- res2$sum
  list(beta[i], mean(res2$sum[501:1000]), sd(res2$sum[501:1000]), matu[,i])
}
t2 <- proc.time()
t2 - t1

save(sample_data, file = "results/potts_sample.rda")
