### Fit the surrogate model

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

d_GP2 <- surrogate_GP(f = mean + var ~ beta - 1, data = df, layers = layers,
                     model = "meanonly",
                     method = "ML", nsteps = 200L,
                     noise_var_m = df$noise_var_m, noise_var_v = df$noise_var_v
)

d_GP3 <- surrogate_GP_stat(f = mean + var ~ beta - 1, data = df,
                           method = "ML", nsteps = 200L,
                           noise_var_m = df$noise_var_m, noise_var_v = df$noise_var_v
)
