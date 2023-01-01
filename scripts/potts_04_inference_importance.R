### Importance sampling with the Potts model

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

##############################

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

for (j in 1:5){
load(paste0('results/potts_simulated_data_', j, '.rda'))

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

log_weights <- summary_stats * beta + 
  summary_stats_sim * beta_hat - summary_stats_sim * beta - log(syn_like)
weights <- exp(log_weights - max(log_weights))
weights_norm <- weights / sum(weights)
posterior_GP <- list(beta, weights_norm, time, time0)

save(posterior_GP, file = paste0('results/potts_importance_', j, '.rda'))
}