### Exchange algorithm with the hidden Potts model

# Load source
source('scripts/utils.R')
load('data/NDVIp089r079_20150503.rda')

# Setup from pixel data
mask <- matrix(1,nrow(ndvi),ncol(ndvi))
n <- prod(dim(mask))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
k <- 5
bcrit <- log(1 + sqrt(k))
y <- as.vector(ndvi)

# Priors
priors <- list()
priors$k <- k
priors$mu <- seq(-0.8, by=0.4, length.out=k)
priors$mu.sd <- rep(0.05,k)
priors$sigma <- rep(0.05,k)
priors$sigma.nu <- rep(k,k)
priors$beta <- c(0.9, 1.3)

### exchange algorithm
## burn-in
nMC <- 2201
nburn <- 202
beta <- rep(0, nMC)
beta[1] <- 1
acc <- 0
sumZs <- numeric(nMC)
muEst <- sdEst <- matrix(nrow=nMC, ncol=priors$k)

# initialize the parameters by sampling from the prior
prSS <- priors$sigma.nu * priors$sigma^2
muEst[1,] <- rnorm(priors$k, priors$mu, priors$mu.sd)
sdEst[1,] <- 1/sqrt(rgamma(priors$k, priors$sigma.nu/2, prSS/2))

# initialize the pixel labels 
n <- length(y)
z <- matrix(0, nrow=n+1, ncol=priors$k)
randZ <- sample.int(priors$k, n, replace = TRUE)
for (i in 1:n) {
  z[i,(randZ[i])] <- 1
}

set.seed(12)

t1 <- proc.time()
for (i in 2:nMC){

  if (i < nburn){bw <- 0.005} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, block, priors, 1)
  sumZs[i] <- summary_stats <- resGibbs$sum
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  ## Simulate
  print(system.time(res2 <- swNoData(beta_star, k, neigh, block, 501)))
  summary_stats_star <- res2$sum[501]
  
  log_MH_ratio <- beta_star * (summary_stats - summary_stats_star) +
    beta[i-1] * (summary_stats_star - summary_stats)
  log_accept_prob <- log(runif(1, 0, 1))
  if (log_MH_ratio > log_accept_prob){
    beta[i] <- beta_star
    if (i >= nburn) {acc <- acc + 1}
  }
  else{
    beta[i] <- beta[i-1]
  }
  print(paste(i, beta[i], sumZs[i], acc))
}
t2 <- proc.time()
time_total <- t2 - t1
time_burnin <- t3 - t1

exchange <- list(beta, muEst, sdEst, sumZs, z, acc, time_total, time_burnin)
save(exchange, file = 'results/hidden_potts_exchange.rda')
