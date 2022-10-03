# Exchange algorithm with the Potts model

# Load source
source('scripts/utils.R')
load("results/potts_sample.rda")

beta <- sample_mean <- sample_var <- rep(0, length(sample_data))

for (i in 1:length(sample_data)){
  beta[i] <- sample_data[[i]][[1]]
  sample_mean[i] <- sample_data[[i]][[2]]
  sample_var[i] <- (sample_data[[i]][[3]])^2
}

df_all <- data.frame(beta = beta, mean = sample_mean, var = sample_var)

df_all$noise_var_m <- df_all$var / 500
df_all$noise_var_v <- 2 * df_all$var^2 / 499

df_test <- df_all[2*(1:50),]
df <- setdiff(df_all, df_test)

mean1 <- mean(df$mean)
df$mean <- df$mean - mean1

# Fit surrogate model

# Gaussian process

layers <- c(AWU(r = 100L, dim = 1L, grad = 250, lims = c(-0.5, 0.5)))

d_GP <- surrogate_GP(f = mean + var ~ beta - 1, data = df, layers = layers,
                     model = "meanvar",
                     method = "ML", nsteps = 200L,
                     noise_var_m = df$noise_var_m, noise_var_v = df$noise_var_v
)

## Load dataset
j = commandArgs(trailingOnly=TRUE)
j <- as.numeric(j)
load(paste0('results/potts_simulated_data_', j, '.rda'))

# Setup
mask <- matrix(1,1e3,1e3)
n <- prod(dim(mask))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
k <- 5
bcrit <- log(1 + sqrt(k))

### exchange algorithm
## burn-in
nMC <- 2201
nburn <- 202
beta <- rep(0, nMC)
beta[1] <- 1
acc <- 0

set.seed(18)

t1 <- proc.time()
for (i in 2:nMC){

  if (i < nburn){bw <- 0.005} else {bw <- 0.001}
  if (i == nburn){t3 <- proc.time()}
  
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
  print(paste(i, beta[i], acc))
}
t2 <- proc.time()
time_total <- t2 - t1
time_burnin <- t3 - t1

exchange <- list(beta, acc, time_total, time_burnin)
save(exchange, file = paste0('results/potts_exchange_', j, '.rda'))
