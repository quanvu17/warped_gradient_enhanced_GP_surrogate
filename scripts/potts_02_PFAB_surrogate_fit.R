### Fit the surrogate model: PFAB algorithm from Moores et al. (2020)

# Load source ----
source('scripts/utils.R')
load("results/potts_sample.rda")

beta <- sample_mean <- sample_var <- rep(0, length(sample_data))
c_idx <- 501:length(sample_data[[1]][[4]])
obsMx <- matrix(nrow=length(sample_data), ncol=length(c_idx))

for (i in 1:length(sample_data)){
  beta[i] <- sample_data[[i]][[1]]
  sample_mean[i] <- sample_data[[i]][[2]]
  sample_var[i] <- (sample_data[[i]][[3]])^2
  obsMx[i,] <- (sample_data[[i]][[4]])[c_idx]
}

df_all <- data.frame(beta = beta, mean = sample_mean, var = sample_var)

df_all$noise_var_m <- df_all$var / 500
df_all$noise_var_v <- 2 * df_all$var^2 / 499

df_test <- df_all[2*(1:50),]
df <- dplyr::setdiff(df_all, df_test)
# plot(df$beta, df$mean, xlab=expression(beta), ylab="S(z)", xlim=c(0,1.3), ylim=c(399000,1998000))

# Initialize parameters ----
n <- 1000
d <- 2
k <- 5

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

bcrit = log(1 + sqrt(k))
print(paste("critical value of beta =",bcrit))
abline(v=bcrit,col=2,lty=3)
E0 <- maxS/k
points(0, E0, col=4, pch=2)

V0 <- maxS*(1/k)*(1 - 1/k)
plot(df$beta, df$var, xlab=expression(beta), ylab=expression(sigma[S(z)]^2),
     pch=20, xlim=c(0.8,1.4),col=4, ylim=c(0,max(df$var, 2*maxS*log(maxS)/pi)))
abline(h=0,lty=2)
abline(v=bcrit,lty=3,col=2)
points(0,V0,col=6,pch=3,lwd=2)
points(bcrit, 2*maxS*log(maxS)/pi, col=6, pch=4, lwd=3)

# Stan ----
library(rstan)
options(mc.cores = min(4,parallel::detectCores()))
M <- length(df$beta)
N <- length(c_idx)
t_idx <- which(beta %in% df$beta)
dat <- list(M=length(df$beta), N=N, maxY=maxS, e0=E0, v0=V0, Vlim=2*maxS*log(maxS)/pi,
            tcrit=bcrit, y=obsMx[t_idx,], t=df$beta)

system.time(fit <- stan("scripts/SurrogatePotts.stan", model_name = "Potts", data = dat, iter=10000,
                        control = list(adapt_delta = 0.99,max_treedepth=20), verbose = TRUE))
# Runtime: 1640.504 seconds
# print(fit, pars=c("a","b","ecrit","vmaxLo","vmaxHi"), digits=3)
save(fit, file=paste0("results/potts_Stan.rda"))
