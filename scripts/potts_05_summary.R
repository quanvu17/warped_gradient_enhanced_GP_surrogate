### Summarize the results of the example with the Potts model

# Load source
source('scripts/utils.R')

set.seed(8)
beta_true <- rnorm(1, log(1 + sqrt(5)), 0.05)

summary_table <- list()

for (j in 1:5){
  ## Surrogate
  load(paste0("results/potts_importance_", j, ".rda"))
  beta_IS <- posterior_GP[[1]]
  # Posterior mean
  post_mean_SR <- mean(beta_IS)
  # Posterior sd
  post_sd_SR <- sd(beta_IS)
  # Time
  time_SR <- posterior_GP[[4]][3] / 3600
  # ESS/time
  ESS_time_SR <- length(beta_IS) / time_SR
  
  ## Importance sampling
  weights_norm <- posterior_GP[[2]]
  # Posterior mean
  post_mean_IS <- sum(weights_norm * beta_IS)
  # Posterior sd
  post_sd_IS <- sqrt(sum(weights_norm * (beta_IS - post_mean_IS)^2))
  # Time
  time_IS <- posterior_GP[[3]][3] / 3600
  # ESS/time
  ESS_time_IS <- sum(weights_norm)^2 / sum(weights_norm^2) / time_IS
  
  ## Delayed-acceptance
  load(paste0("results/potts_delayed_accept_", j, ".rda"))
  beta_DA <- delay_accept[[1]][202:2201]
  # Posterior mean
  post_mean_DA <- mean(beta_DA)
  # Posterior sd
  post_sd_DA <- sd(beta_DA)
  # Time
  time_DA <- (delay_accept[[4]][3] - delay_accept[[5]][3]) / 3600
  # ESS/time
  ESS_time_DA <- coda::effectiveSize(beta_DA) / time_DA
  
  ## Exchange
  load(paste0("results/potts_exchange_", j, ".rda"))
  beta_EX <- exchange[[1]][202:2201]
  # Posterior mean
  post_mean_EX <- mean(beta_EX)
  # Posterior sd
  post_sd_EX <- sd(beta_EX)
  # Time
  time_EX <- (exchange[[3]][3] - exchange[[4]][3]) / 3600
  # ESS/time
  ESS_time_EX <- coda::effectiveSize(beta_EX) / time_EX
  
  summary_table[[j]] <- matrix(c(post_mean_SR, post_sd_SR, time_SR, ESS_time_SR,
                                 post_mean_IS, post_sd_IS, time_IS, ESS_time_IS,
                                 post_mean_DA, post_sd_DA, time_DA, ESS_time_DA,
                                 post_mean_EX, post_sd_EX, time_EX, ESS_time_EX),
                               nrow = 4)
  
  ## Plot
  
  df_IS <- data.frame(beta = beta_IS, weights_norm = weights_norm)
  df_DA <- data.frame(beta = beta_DA)
  df_EX <- data.frame(beta = beta_EX)
  
  plot <- ggplot() + 
    stat_density(data = df_IS, aes(beta), adjust = 1.5, color = "purple", geom = "line", linetype = "dashed", size = 0.8) +
    stat_density(data = df_IS, aes(beta, weight = weights_norm), adjust = 1.5, color = "blue", geom = "line", linetype = "dotted", size = 0.8) +
    stat_density(data = df_DA, aes(beta), adjust = 1.5, color = "darkgreen", geom = "line", linetype = "dotdash", size = 0.8) +
    stat_density(data = df_EX, aes(beta), adjust = 1.5, color = "red", geom = "line", size = 0.8) +
    geom_vline(xintercept = beta_true, colour = "black", size = 0.5, linetype="dashed") +
    theme_bw() + labs(x=expression(beta^'(1)'), y='density')
  
  ggsave(filename=paste0("figures/potts_posterior_density_plot_", j, ".png"),
         plot=plot, device="png", width=10, height=10, scale=1, units="cm", dpi=300)
  
  
}