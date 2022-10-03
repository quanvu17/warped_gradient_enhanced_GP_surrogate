### Summarize the results of the example with the autologistic model

# Load source
source('scripts/utils.R')

# Importance sampling
load("results/autologistic_importance.rda")
beta_IS <- posterior_GP[[1]]
# Posterior mean
post_mean1_SR <- mean(beta_IS[,1])
post_mean2_SR <- mean(beta_IS[,2])
# Posterior sd
post_sd1_SR <- sd(beta_IS[,1])
post_sd2_SR <- sd(beta_IS[,2])
# Time
time_SR <- posterior_GP[[4]][3] / 60
# ESS/time
ESS_time_SR <- nrow(beta_IS) / time_SR

# Importance sampling
weights_norm <- posterior_GP[[2]]
# Posterior mean
post_mean1_IS <- sum(beta_IS[,1] * weights_norm)
post_mean2_IS <- sum(beta_IS[,2] * weights_norm)
# Posterior sd
post_sd1_IS <- sqrt(sum(weights_norm * (beta_IS[,1] - post_mean1_IS)^2))
post_sd2_IS <- sqrt(sum(weights_norm * (beta_IS[,2] - post_mean2_IS)^2))
# Time
time_IS <- posterior_GP[[3]][3] / 60
# ESS/time
ESS_time_IS <- sum(weights_norm)^2 / sum(weights_norm^2) / time_IS

# Exchange
load("results/autologistic_exchange.rda")
beta_EX <- exchange[[1]][202:4201,]
# Posterior mean
post_mean1_EX <- mean(beta_EX[,1])
post_mean2_EX <- mean(beta_EX[,2])
# Posterior sd
post_sd1_EX <- sd(beta_EX[,1])
post_sd2_EX <- sd(beta_EX[,2])
# Time
time_EX <- (exchange[[3]] - exchange[[4]])[3] / 60
# ESS/time
ESS_time_EX <- coda::effectiveSize(beta_EX) / time_EX

### Plot
df_IS <- data.frame(beta1 = beta_IS[,1], beta2 = beta_IS[,2], weights_norm = weights_norm)
df_EX <- data.frame(beta1 = beta_EX[,1], beta2 = beta_EX[,2])

plot1 <- ggplot() + 
  stat_density(data = df_IS, aes(beta1, weight = weights_norm), adjust = 2, color = "blue", geom = "line", linetype = "dotted", size = 0.8) +
  stat_density(data = df_IS, aes(beta1), adjust = 2, color = "purple", geom = "line", linetype = "dashed", size = 0.8) +
  stat_density(data = df_EX, aes(beta1), adjust = 2, color = "red", geom = "line", size = 0.8) +
  xlim(-0.015, 0.01) +
  theme_bw() + labs(x=expression(beta^'(1)'), y='density')

plot2 <- ggplot() + 
  stat_density(data = df_IS, aes(beta2, weight = weights_norm), adjust = 2, color = "blue", geom = "line", linetype = "dotted", size = 0.8) +
  stat_density(data = df_IS, aes(beta2), adjust = 2, color = "purple", geom = "line", linetype = "dashed", size = 0.8) +
  stat_density(data = df_EX, aes(beta2), adjust = 2, color = "red", geom = "line", size = 0.8) +
  xlim(0.80, 0.95) +
  theme_bw() + labs(x=expression(beta^'(2)'), y='density')

plot <- grid.arrange(plot1, plot2, nrow = 1)

ggsave(filename="figures/autologistic_posterior_density_plot.png", plot=plot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)

