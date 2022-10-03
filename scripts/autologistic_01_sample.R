### Simulate the pseudo-data to fit the surrogate model

# Load source
source('scripts/utils.R')

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

registerDoParallel(cores = 64)
nsim = 30

## Step 1
# Optional potential on sites
beta1 <- seq(-0.2, 0.1, by = 0.05) # seq(-0.4, 0.1, by = 0.05)
# Interaction parameter
beta2 <- seq(0.7, 1.2, by = 0.05)
beta <- expand.grid(beta1 = beta1, beta2 = beta2)

# Sampling using an existing configuration as starting point
mean1 <- mean2 <- var1 <- var2 <- list()
autologistic <- foreach (j = 1:nrow(beta)) %dopar% {
  summarystats1 <- summarystats2 <- list()
  for (i in 1:nsim){
    set.seed(i)
    image <- sampler.mrf(iter = n, sampler = method, h = height, w = width,
                         ncolors = K, nei = G, param = beta[j,2],
                         pot = c(beta[j,1], -beta[j,1]), ## Ising labels 1 and -1
                         initialise = FALSE, view = FALSE)
    summarystats1[[i]] <- sum(image == 0) - sum(image == 1)
    summarystats2[[i]] <- CountNeighbors(image)
  }
  mean1 <- mean(unlist(summarystats1))
  var1 <- var(unlist(summarystats1))
  mean2 <- mean(unlist(summarystats2))
  var2 <- var(unlist(summarystats2))
  c(mean1, mean2, var1, var2)
}

df <- matrix(unlist(autologistic), ncol=4, byrow=T)
df <- data.frame(beta[,1], beta[,2], df)
names(df)[1:6] <- c("beta1", "beta2", "mean_1", "mean_2", "var_1", "var_2")
save(df, file = "results/autologistic_sample.rda")
