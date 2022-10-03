### Summarize the results of the example with the hidden Potts model

# Load source
source('scripts/utils.R')

## Surrogate
load("results/hidden_potts_surrogate.rda")
beta_SR <- surrogate[[1]][202:2201]
# Posterior mean
post_mean_SR <- mean(beta_SR)
# Posterior sd
post_sd_SR <- sd(beta_SR)
# Time
time_SR <- (surrogate[[6]][3] - surrogate[[7]][3]) / 3600
# ESS/time
ESS_time_SR <- coda::effectiveSize(beta_SR) / time_SR

## Delayed-acceptance
load("results/hidden_potts_delayed_accept.rda")
beta_DA <- delay_accept[[1]][202:2201]
# Posterior mean
post_mean_DA <- mean(beta_DA)
# Posterior sd
post_sd_DA <- sd(beta_DA)
# Time
time_DA <- (delay_accept[[8]][3] - delay_accept[[9]][3]) / 3600
# ESS/time
ESS_time_DA <- coda::effectiveSize(beta_DA) / time_DA

## Exchange
load("results/hidden_potts_exchange.rda")
beta_EX <- exchange[[1]][202:2201]
# Posterior mean
post_mean_EX <- mean(beta_EX)
# Posterior sd
post_sd_EX <- sd(beta_EX)
# Time
time_EX <- (exchange[[7]][3] - exchange[[8]][3]) / 3600
# ESS/time
ESS_time_EX <- coda::effectiveSize(beta_EX) / time_EX

## Plot

df_SR <- data.frame(beta = beta_SR)
df_DA <- data.frame(beta = beta_DA)
df_EX <- data.frame(beta = beta_EX)

plot <- ggplot() + 
  stat_density(data = df_SR, aes(beta), adjust = 2, color = "purple", geom = "line", linetype = "dashed", size = 0.8) +
  stat_density(data = df_DA, aes(beta), adjust = 2, color = "darkgreen", geom = "line", linetype = "dotdash", size = 0.8) +
  stat_density(data = df_EX, aes(beta), adjust = 2, color = "red", geom = "line", size = 0.8) +
  xlim(c(1.232, 1.239)) +
  theme_bw() + labs(x=expression(beta^'(1)'), y='density')

ggsave(filename="figures/hidden_potts_posterior_density_plot.png", plot=plot, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

