### Exchange algorithm with the autologistic model

# Load source
source('scripts/utils.R')

## Summary stats for the IceFloe image
data("icefloe")
summarystats1 <- sum(IceFloe == F) - sum(IceFloe == T)
summarystats2 <- CountNeighbors(IceFloe)

# Algorithm settings
n <- 200
method <- "SW"
# Dimension of the lattice
height <- width <- 40
# Number of colors. Automatically set to 2 if not specified.
K <- 2
# Number of neighbors. Automatically set to 4 if not specified.
G <- 4

### exchange algorithm
nMC <- 4201
nburn <- 202
beta <- matrix(rep(0, 2*nMC), ncol = 2)
beta[1,1] <- -0.05
beta[1,2] <- 0.7
acc <- 0

set.seed(27)

t1 <- proc.time()
for (i in 2:nMC){
  
  if (i < nburn){
    bw1 <- 0.02
    bw2 <- 0.1
  } else {
    bw1 <- 0.005
    bw2 <- 0.025
  }
  if (i == nburn){t3 <- proc.time()}
  
  beta1_star <- rnorm(1, beta[i-1, 1], bw1)
  beta2_star <- rnorm(1, beta[i-1, 2], bw2)
  
  ## Simulate
  image <- sampler.mrf(iter = n, sampler = method, h = height, w = width,
                       ncolors = K, nei = G, param = beta2_star,
                       pot = c(beta1_star, -beta1_star), ## Ising labels 1 and -1
                       initialise = FALSE, view = FALSE)
  summarystats1_star <- sum(image == 0) - sum(image == 1)
  summarystats2_star <- CountNeighbors(image)
  
  log_MH_ratio <- (beta1_star - beta[i-1, 1]) * (summarystats1 - summarystats1_star) +
    (beta2_star - beta[i-1, 2]) * (summarystats2 - summarystats2_star)
    
  log_accept_prob <- log(runif(1, 0, 1))
  
  if (log_MH_ratio > log_accept_prob){
    beta[i, 1] <- beta1_star
    beta[i, 2] <- beta2_star
    if (i >= nburn) {acc <- acc + 1}
  }
  else{
    beta[i, 1] <- beta[i-1, 1]
    beta[i, 2] <- beta[i-1, 2]
  }
  print(paste(i, beta[i, 1], beta[i, 2], acc))
}
t2 <- proc.time()
time_total <- t2 - t1
time_burnin <- t3 - t1

exchange <- list(beta, acc, time_total, time_burnin)
save(exchange, file = 'results/autologistic_exchange.rda')
