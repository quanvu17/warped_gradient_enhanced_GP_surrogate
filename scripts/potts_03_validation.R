### Compare the different surrogate models

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

newdata <- data.frame(beta = df_test$beta)
load('results/potts_simulated_data.rda')

syn_like_GP <- synthetic_likefnc_GP_1D(d_GP, newdata, summary_stats, mean1)
syn_like_GP2 <- synthetic_likefnc_GP_1D(d_GP2, newdata, summary_stats, mean1)
syn_like_GP3 <- synthetic_likefnc_GP_1D_stat(d_GP3, newdata, summary_stats, mean1)

## Calculate MAE
MAE1 <- mean(abs(df_test$mean - syn_like_GP[[1]]))
MAE2 <- mean(abs(df_test$mean - syn_like_GP2[[1]]))
MAE3 <- mean(abs(df_test$mean - syn_like_GP3[[1]]))

## Calculate RMSPE
RMSE1 <- sqrt(mean((df_test$mean - syn_like_GP[[1]])^2))
RMSE2 <- sqrt(mean((df_test$mean - syn_like_GP2[[1]])^2))
RMSE3 <- sqrt(mean((df_test$mean - syn_like_GP3[[1]])^2))

## Calculate MAE
MAE1_v <- mean(abs(df_test$var - syn_like_GP[[2]]))
MAE2_v <- mean(abs(df_test$var - syn_like_GP2[[2]]))
MAE3_v <- mean(abs(df_test$var - syn_like_GP3[[2]]))

## Calculate RMSPE
RMSE1_v <- sqrt(mean((df_test$var - syn_like_GP[[2]])^2))
RMSE2_v <- sqrt(mean((df_test$var - syn_like_GP2[[2]])^2))
RMSE3_v <- sqrt(mean((df_test$var - syn_like_GP3[[2]])^2))

pred <- c(MAE1, MAE2, MAE3, RMSE1, RMSE2, RMSE3,
          MAE1_v, MAE2_v, MAE3_v, RMSE1_v, RMSE2_v, RMSE3_v)
save(pred, file = "results/potts_prediction_results.rda")

plot <- ggplot(data = df_test) + 
  geom_point(aes(beta, abs(mean - syn_like_GP[[1]])), col = "blue") +
  geom_point(aes(beta, abs(mean - syn_like_GP2[[1]])), col = "red", shape = 15) + 
  geom_point(aes(beta, abs(mean - syn_like_GP3[[1]])), col = "darkgreen", shape = 17) +
  xlim(c(1.05, 1.25)) +
  theme_bw() + labs(x=expression(beta^'(1)'), y='absolute error')
  
ggsave(filename="figures/potts_absolute_error_plot.png", plot=plot, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

## Plot

plot1 <- ggplot() + 
  geom_point(data = df, aes(beta, mean + mean1)) +
  geom_point(data = df_test, aes(beta, mean)) + # shape = 2
  xlim(c(1.05, 1.25)) +
  theme_bw() + labs(x=expression(beta^'(1)'), y='mean')

plot2 <- ggplot() + 
  geom_point(data = df, aes(beta, sqrt(var))) +
  geom_point(data = df_test, aes(beta, sqrt(var))) + # shape = 2
  xlim(c(1.05, 1.25)) +
  theme_bw() + labs(x=expression(beta^'(1)'), y='standard deviation')

plot <- grid.arrange(plot1, plot2, nrow = 1)

ggsave(filename="figures/potts_summary_stats_plot.png", plot=plot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)
