### Delayed-acceptance MCMC with the Potts model

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

# Load source
i = commandArgs(trailingOnly=TRUE)
i <- as.numeric(i)

load(paste0('results/potts_simulated_data_', i, '.rda'))

# Setup
mask <- matrix(1,1e3,1e3)
n <- prod(dim(mask))
neigh <- getNeighbors(mask, c(2,2,0,0))
block <- getBlocks(mask, 2)
edges <- getEdges(mask, c(2,2,0,0))
maxS <- nrow(edges)
k <- 5
bcrit <- log(1 + sqrt(k))

### delayed acceptance algorithm
nMC <- 2201
nburn <- 202
beta <- syn_like <- rep(0, nMC)
beta[1] <- 1

newdata <- data.frame(beta = beta[1])
syn_like_GP <- synthetic_likefnc_GP_1D(d_GP, newdata, summary_stats, mean1)
syn_like[[1]] <- syn_like_GP[[3]]

acc <- 0
acc_step1 <- 0

t1 <- proc.time()
for (i in 2:nMC){
  
  if (i < nburn){bw <- 0.05} else {bw <- 0.001}
  if (i == nburn){t3 <- proc.time()}
  
  set.seed(i)
  beta_star <- rnorm(1, beta[i-1], bw)
  
  # Step 1: use the surrogate likelihood
  if (i < nburn){
    newdata <- data.frame(beta = beta_star)
    syn_like_GP <- synthetic_likefnc_GP_1D(d_GP, newdata, summary_stats, mean1)
    syn_like_star <- syn_like_GP[[3]]
  } else {
    newdata <- data.frame(beta = beta_star)
    syn_like_GP <- synthetic_likefnc_GP_1D_sample(d_GP, newdata, summary_stats, mean1)
    syn_like_star <- log(rowMeans(exp(syn_like_GP[[4]]), na.rm = T))
  }
  
  log_MH_ratio <- syn_like_star - syn_like[i-1]
  set.seed(i)
  log_accept_prob <- log(runif(1, 0, 1))
  
  if (i < nburn){
    if (log_MH_ratio > log_accept_prob){
      beta[i] <- beta_star
      syn_like[i] <- syn_like_star
    }
    else{
      beta[i] <- beta[i-1]
      syn_like[i] <- syn_like[i-1]
    }
  } else {
    if (log_MH_ratio > log_accept_prob){ ## Proceed to Step 2
      
      acc_step1 <- acc_step1 + 1
      
      ## Simulate
      print(system.time(res2 <- swNoData(beta_star, k, neigh, block, 501)))
      summary_stats_star <- res2$sum[501]
      
      log_MH_ratio_2 <- - syn_like_star + syn_like[i-1] + 
        beta_star * (summary_stats - summary_stats_star) +
        beta[i-1] * (summary_stats_star - summary_stats)
      log_accept_prob_2 <- log(runif(1, 0, 1))
      if (log_MH_ratio_2 > log_accept_prob_2){
        beta[i] <- beta_star
        syn_like[i] <- syn_like_star
        acc <- acc + 1
      }
      else{
        beta[i] <- beta[i-1]
        syn_like[i] <- syn_like[i-1]
      }
    } 
    else{
      beta[i] <- beta[i-1]
      syn_like[i] <- syn_like[i-1]
    }
  }
  print(paste(i, beta[i], acc_step1, acc))
}
t2 <- proc.time()
time_total <- t2 - t1
time_burnin <- t3 - t1

delay_accept <- list(beta, acc, acc_step1, time_total, time_burnin)
save(delay_accept, file = paste0('results/potts_delayed_accept_', i, '.rda'))
