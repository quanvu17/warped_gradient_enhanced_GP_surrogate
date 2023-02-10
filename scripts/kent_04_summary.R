### Summarize the results of the example with the Kent model

# Load source
source('scripts/utils.R')
load("results/kent_importance.rda")

param_IS <- posterior_GP[[1]]
# Importance sampling
weights_norm <- posterior_GP[[2]]
# Posterior mean
post_mean1_IS <- sum(param_IS[,1] * weights_norm)
post_mean2_IS <- sum(param_IS[,2] * weights_norm)
# Posterior sd
post_sd1_IS <- sqrt(sum(weights_norm * (param_IS[,1] - post_mean1_IS)^2))
post_sd2_IS <- sqrt(sum(weights_norm * (param_IS[,2] - post_mean2_IS)^2))
# Time
time_IS <- posterior_GP[[3]][3] / 60
# ESS/time
ESS_time_IS <- sum(weights_norm)^2 / sum(weights_norm^2) / time_IS


### Plot
df_IS <- data.frame(param1 = param_IS[,1], param2 = param_IS[,2], weights_norm = weights_norm)

plot1 <- ggplot() + 
  stat_density(data = df_IS, aes(param1, weight = weights_norm), adjust = 2, color = "blue", geom = "line", linetype = "dotted", size = 0.8) +
  stat_density(data = kappa_posterior_unnorm, aes(kappa, weight = dens), adjust = 0.5, color = "black", geom = "line", size = 0.5) +
  geom_vline(xintercept = 2, colour = "black", size = 0.5, linetype="dashed") +
  theme_bw() + labs(x=expression(kappa), y='density') + xlim(1, 3)

plot2 <- ggplot() + 
  stat_density(data = df_IS, aes(param2, weight = weights_norm), adjust = 2, color = "blue", geom = "line", linetype = "dotted", size = 0.8) +
  stat_density(data = beta_posterior_unnorm, aes(beta, weight = dens), adjust = 0.5, color = "black", geom = "line", size = 0.5) +
  geom_vline(xintercept = 0.5, colour = "black", size = 0.5, linetype="dashed") +
  theme_bw() + labs(x=expression(beta), y='density') + xlim(0, 1.5)

plot <- grid.arrange(plot1, plot2, nrow = 1)

ggsave(filename="figures/kent_posterior_density_plot.png", plot=plot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)

