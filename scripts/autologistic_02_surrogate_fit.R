### Fit the surrogate model

# Load source
source('scripts/utils.R')
load("results/autologistic_sample.rda")

mean1 <- mean(df$mean_1)
mean2 <- mean(df$mean_2)

df$mean_1 <- df$mean_1 - mean1
df$mean_2 <- df$mean_2 - mean2

df$noise_var_m_1 <- df$var_1 / 30
df$noise_var_m_2 <- df$var_2 / 30
df$noise_var_v_1 <- 2 * df$var_1^2 / 29
df$noise_var_v_2 <- 2 * df$var_2^2 / 29

### GP

layers <- c(AWU(r = 100L, dim = 1L, grad = 250, lims = c(-0.6, 0.6)),
            AWU(r = 100L, dim = 2L, grad = 250, lims = c(-0.6, 0.6))
)


d_GP_1 <- surrogate_GP_2D(f = mean_1 + var_1 ~ beta1 + beta2 - 1, data = df, layers = layers,
                       model = "meanvar",
                       method = "ML", nsteps = 150L, # 200L
                       noise_var_m = df$noise_var_m_1, noise_var_v = df$noise_var_v_1
)


layers <- c(AWU(r = 100L, dim = 1L, grad = 250, lims = c(-0.6, 0.6)),
            AWU(r = 100L, dim = 2L, grad = 250, lims = c(-0.6, 0.6))
)

d_GP_2 <- surrogate_GP_2D(f = mean_2 + var_2 ~ beta2 + beta1 - 1, data = df, layers = layers,
                       model = "meanvar",
                       method = "ML", nsteps = 150L, # 200L
                       noise_var_m = df$noise_var_m_2, noise_var_v = df$noise_var_v_2
)
