### MCMC (using the surrogate likelihood) with the hidden Potts model

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

### delayed acceptance algorithm
nMC <- 2201
nburn <- 202
beta <- syn_like <- rep(0, nMC)
beta[1] <- 1

newdata <- data.frame(beta = beta[1])

acc <- 0
acc_step1 <- 0

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

set.seed(23)

t1 <- proc.time()
for (i in 2:nMC){

  if (i < nburn){bw <- 0.05} else {bw <- 0.002}
  if (i == nburn){t3 <- proc.time()}
  
  # Gibbs sampler for the mu, sigma parameters, and update summary stats
  resGibbs <- gibbsPotts(y, z, beta[i-1], muEst[i-1,], sdEst[i-1,], neigh, block, priors, 1)
  sumZs[i] <- summary_stats <- as.numeric(resGibbs$sum)
  muEst[i,] <- resGibbs$mu
  sdEst[i,] <- resGibbs$sigma
  z <- resGibbs$z
  
  beta_star <- rnorm(1, beta[i-1], bw)
  
  # Step 1: use the surrogate likelihood
  newdata <- data.frame(beta = c(beta[i-1], beta_star))
  syn_like_GP <- synthetic_likefnc_GP_1D(d_GP, newdata, summary_stats, mean1)
  syn_like[i-1] <- syn_like_GP[[3]][1]
  syn_like_star <- syn_like_GP[[3]][2]
  
  log_MH_ratio <- syn_like_star - syn_like[i-1]
  log_accept_prob <- log(runif(1, 0, 1))
  
    if (log_MH_ratio > log_accept_prob){
      beta[i] <- beta_star
    }
    else{
      beta[i] <- beta[i-1]
    }
  
  print(paste(i, beta[i], sumZs[i]))
}
t2 <- proc.time()
time_total <- t2 - t1
time_burnin <- t3 - t1

surrogate <- list(beta, muEst, sdEst, sumZs, z, time_total, time_burnin)
save(surrogate, file = 'results/hidden_potts_surrogate.rda')
