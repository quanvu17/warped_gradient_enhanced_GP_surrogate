### Importance sampling (using the surrogate model) with the Kent model

# Load source
source('scripts/utils.R')
library(Directional)

registerDoParallel(cores = detectCores()) ## 64L

## Summary stats
load('results/kent_simulated_data.rda')
load("results/kent_sample.rda")
mean1 <- mean(df$mean_1)
mean2 <- mean(df$mean_2)

############## GP

k <- seq(1.2, 3, by = 0.025)
b <- seq(0.1, 1.5, by = 0.0125)
param <- expand.grid(k = k, b = b)
param <- param %>% filter(b/k < 0.5)
newdata <- data.frame(param1 = param[,1], param2 = param[,2])
newdata.1 <- newdata[1:2774,]
newdata.2 <- newdata[2775:5548,]

syn_like_GP.1 <- synthetic_likefnc_GP_2D(d_GP_1, d_GP_2, newdata.1, summary_stats, mean1, mean2)
syn_like_GP.2 <- synthetic_likefnc_GP_2D(d_GP_1, d_GP_2, newdata.2, summary_stats, mean1, mean2)

syn_like.1 <- exp(syn_like_GP.1[[5]])
syn_like.2 <- exp(syn_like_GP.2[[5]])

syn_like <- c(syn_like.1, syn_like.2)

### Importance Sampling
param_hat <- c(newdata$param1[which.max(syn_like)],
               newdata$param2[which.max(syn_like)]) ## estimated maximum likelihood estimator

set.seed(2022)
idx <- sample(1:nrow(newdata), size = 2000, prob = syn_like / sum(syn_like), replace = T)
syn_like <- syn_like[idx]
param <- matrix(c(newdata$param1[idx],
                  newdata$param2[idx]), byrow = F, ncol = 2)
plot(param)

## simulate
t1 <- proc.time()
mean1 <- mean2 <- var1 <- var2 <- list()
kent <- foreach (j = 1:nrow(param)) %dopar% {
  summarystats1 <- summarystats2 <- list()
  y <- rkent2(1e2, param[j,1], gamma_1, param[j,2])
  summarystats1 <- sum(y[,1] * gamma_1[1])
  summarystats2 <- sum((y[,2] * gamma_2[2])^2 - (y[,3] * gamma_3[3])^2)
  c(summarystats1, summarystats2)
}
t2 <- proc.time()
time <- t2 - t1

summary_stats_sim1 <- summary_stats_sim2 <- rep(0, nrow(param))

for (i in 1:nrow(param)){
  summary_stats_sim1[i] <- kent[[i]][[1]]
  summary_stats_sim2[i] <- kent[[i]][[2]]
}

log_weights <- summary_stats[1] * param[,1] + summary_stats[2] * param[,2] +
  summary_stats_sim1 * param_hat[1] + summary_stats_sim2 * param_hat[2] +
  (- summary_stats_sim1 * param[,1]) +  (- summary_stats_sim2 * param[,2]) +
  (- log(syn_like))

weights <- exp(log_weights - max(log_weights))
weights_norm <- weights / sum(weights)
posterior_GP <- list(param, weights_norm, time)

save(posterior_GP, file = "results/kent_importance.rda")
