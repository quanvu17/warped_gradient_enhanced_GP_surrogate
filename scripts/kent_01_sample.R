### Simulate the pseudo-data to fit the surrogate model

# Load source
source('scripts/utils.R')
library(Directional)

registerDoParallel(cores = 64)
nsim = 100

## simulate from the Kent distribution

gamma_1 <- c(1, 0, 0)
gamma_2 <- c(0, 1, 0)
gamma_3 <- c(0, 0, 1)

# Parameter
k <- seq(0.2, 10, by = 0.2)
b <- seq(0.1, 5, by = 0.1)
param <- expand.grid(k = k, b = b)
param <- param %>% filter(b/k < 0.5)

# Sampling
mean1 <- mean2 <- var1 <- var2 <- list()
kent <- foreach (j = 1:nrow(param)) %dopar% {
  summarystats1 <- summarystats2 <- list()
  for (i in 1:nsim){
    set.seed(i)
    y <- rkent2(1e2, param$k[j], gamma_1, param$b[j])
    summarystats1[[i]] <- sum(y[,1] * gamma_1[1])
    summarystats2[[i]] <- sum((y[,2] * gamma_2[2])^2 - (y[,3] * gamma_3[3])^2)
    # summarystats1[[i]] <- sum(y[,1] * gamma_1[1] + y[,2] * gamma_1[2] + y[,3] * gamma_1[3])
    # summarystats2[[i]] <- sum((y[,1] * gamma_2[1] + y[,2] * gamma_2[2] + y[,3] * gamma_2[3])^2 -
    #                          (y[,1] * gamma_3[1] + y[,2] * gamma_3[2] + y[,3] * gamma_3[3])^2)
  }
  mean1 <- mean(unlist(summarystats1))
  var1 <- var(unlist(summarystats1))
  mean2 <- mean(unlist(summarystats2))
  var2 <- var(unlist(summarystats2))
  c(mean1, mean2, var1, var2)
}

df <- matrix(unlist(kent), ncol=4, byrow=T)
df <- data.frame(param[,1], param[,2], df)
names(df)[1:6] <- c("param1", "param2",
                    "mean_1", "mean_2", "var_1", "var_2")
save(df, file = "results/kent_sample.rda")
