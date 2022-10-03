### Importance sampling with the Potts model

# Load source
source('scripts/utils.R')
registerDoParallel(cores = detectCores()) ## 64L

# Initialize parameters
n <- 1000
d <- 2
k <- 5
iter <- 1000

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

############## GP
newdata <- data.frame(beta = seq(1.16, 1.18, length.out = 1001))

i = commandArgs(trailingOnly=TRUE)
i <- as.numeric(i)
load(paste0('results/potts_simulated_data_', i, '.rda'))

syn_like_GP <- synthetic_likefnc_GP_1D_sample(d_GP, newdata, summary_stats, mean1)
syn_like <- rowMeans(exp(syn_like_GP[[4]]), na.rm = T)

### Importance Sampling
beta_hat <- newdata$beta[which.max(syn_like)] ## estimated maximum likelihood estimator

set.seed(2)

t1 <- proc.time()
idx <- sample(1:nrow(newdata), size = 1000, prob = syn_like / sum(syn_like),
              replace = T)
beta <- newdata$beta[idx] + runif(1000, -1e-5, 1e-5)
t2 <- proc.time()
time0 <- t2 - t1

syn_like <- syn_like[idx]

t1 <- proc.time()
sample_data <- foreach (i = 1:length(beta)) %dopar% {
  res2 <- swNoData(beta[i], k, neigh, block, 501)
  list(beta[i], res2$sum[501])
}
t2 <- proc.time()
time <- t2 - t1

summary_stats_sim <- rep(0, length(beta))

for (i in 1:length(beta)){
  summary_stats_sim[i] <- sample_data[[i]][[2]]
}

log_weights <- summary_stats * beta + summary_stats_sim * beta_hat - summary_stats_sim * beta - syn_like
weights <- exp(log_weights - max(log_weights))
weights_norm <- weights / sum(weights)
posterior_GP <- list(beta, weights_norm, time, time0)

save(posterior_GP, file = paste0('results/potts_importance_', i, '.rda'))
