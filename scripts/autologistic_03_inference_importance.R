### Importance sampling (using the surrogate model) with the autologistic model

# Load source
source('scripts/utils.R')

registerDoParallel(cores = detectCores()) ## 64L

## Summary stats for the IceFloe image
load("data/data_icefloe.rda")
summarystats1 <- sum(IceFloe == F) - sum(IceFloe == T)
summarystats2 <- CountNeighbors(IceFloe)
summary_stats <- c(summarystats1, summarystats2)

############## GP

beta1 <- seq(-0.02, 0.01, by = 0.001) 
beta2 <- seq(0.8, 0.95, by = 0.001)
beta <- expand.grid(beta1 = beta1, beta2 = beta2)
newdata <- data.frame(beta1 = beta[,1], beta2 = beta[,2])
newdata.1 <- newdata[1:2300,]
newdata.2 <- newdata[2301:4681,]

syn_like_GP.1 <- synthetic_likefnc_GP_2D(d_GP_1, d_GP_2, newdata.1, summary_stats, mean1, mean2)
syn_like_GP.2 <- synthetic_likefnc_GP_2D(d_GP_1, d_GP_2, newdata.2, summary_stats, mean1, mean2)

syn_like.1 <- exp(syn_like_GP.1[[5]])
syn_like.2 <- exp(syn_like_GP.2[[5]])

syn_like <- c(syn_like.1, syn_like.2)

### Importance Sampling
beta_hat <- c(newdata$beta1[which.max(syn_like)],
              newdata$beta2[which.max(syn_like)]) ## estimated maximum likelihood estimator

set.seed(2)
t1 <- proc.time()
idx <- sample(1:nrow(newdata), size = 2000, prob = syn_like / sum(syn_like),
              replace = T)
beta <- matrix(c(newdata$beta1[idx] + runif(2000, -0.001/2, 0.001/2),
                 newdata$beta2[idx] + runif(2000, -0.001/2, 0.001/2)), byrow = F, ncol = 2)
t2 <- proc.time()
time0 <- t2 - t1

syn_like <- syn_like[idx]

# sampler.mrf
# Algorithm settings
n <- 200
method <- "SW"
# Dimension of the lattice
height <- width <- 40
# Number of colors. Automatically set to 2 if not specified.
K <- 2
# Number of neighbors. Automatically set to 4 if not specified.
G <- 4

t1 <- proc.time()
mean1 <- mean2 <- var1 <- var2 <- list()
autologistic <- foreach (j = 1:nrow(beta)) %dopar% {
  summarystats1 <- summarystats2 <- list()
    image <- GiRaF::sampler.mrf(iter = n, sampler = method, h = height, w = width,
                         ncolors = K, nei = G, param = beta[j,2],
                         pot = c(beta[j,1], -beta[j,1]), ## Ising labels 1 and -1
                         initialise = FALSE, view = TRUE)
    summarystats1 <- sum(image == 0) - sum(image == 1)
    summarystats2 <- CountNeighbors(image)
  c(summarystats1, summarystats2)
}
t2 <- proc.time()
time <- t2 - t1

summary_stats_sim1 <- summary_stats_sim2 <- rep(0, nrow(beta))

for (i in 1:nrow(beta)){
  summary_stats_sim1[i] <- autologistic[[i]][[1]]
  summary_stats_sim2[i] <- autologistic[[i]][[2]]
}

log_weights <- summary_stats[1] * beta[,1] + summary_stats[2] * beta[,2] +
  summary_stats_sim1 * beta_hat[1] + summary_stats_sim2 * beta_hat[2] +
  (- summary_stats_sim1 * beta[,1]) +  (- summary_stats_sim2 * beta[,2]) +
  (- log(syn_like))

weights <- exp(log_weights - max(log_weights))
weights_norm <- weights / sum(weights)
posterior_GP <- list(beta, weights_norm, time, time0)

save(posterior_GP, file = "results/autologistic_importance.rda")
