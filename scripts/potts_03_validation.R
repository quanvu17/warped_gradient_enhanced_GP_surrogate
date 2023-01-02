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
load('results/potts_simulated_data_1.rda')

# Prediction for GP surrogate models
syn_like_GP <- synthetic_likefnc_GP_1D(d_GP, newdata, summary_stats, mean1)
syn_like_GP2 <- synthetic_likefnc_GP_1D(d_GP2, newdata, summary_stats, mean1)
syn_like_GP3 <- synthetic_likefnc_GP_1D_stat(d_GP3, newdata, summary_stats, mean1)

# Prediction for PFAB surrogate model
load("results/potts_Stan.rda")
# Predict mean
ft <- function(t, tC, e0, ecrit, v0, vmax1, vmax2, phi1, phi2) {
  sqrtBcritPhi = sqrt(tC)*phi1
  fval <- numeric(length(t))
  for (i in 1:length(t)) {
    if (t[i] <= tC) {
      sqrtBdiffPhi = sqrt(tC - t[i])*phi1
      fval[i] <- e0 + t[i]*v0 - ((2*(vmax1-v0))/(phi1^2))*((sqrtBcritPhi + 1)/exp(sqrtBcritPhi) - (sqrtBdiffPhi + 1)/exp(sqrtBdiffPhi))
    } else {
      sqrtBdiff = sqrt(t[i] - tC)
      fval[i] <- ecrit - ((2*vmax2)/phi2)*(sqrtBdiff/exp(phi2*sqrtBdiff) + (exp(-phi2*sqrtBdiff) - 1)/phi2);
    }
  }
  return(fval)
}
# Predict variance
dfdt <- function(t, tC, V0, Vmax1, Vmax2, r1, r2) {
  ifelse(t < tC,
         V0 + (Vmax1-V0)*exp(-r1*sqrt(tC - t)),
         Vmax2*exp(-r2*sqrt(t - tC)))
}
# Prediction
para <- get_posterior_mean(fit, pars=c("a","b","ecrit","vmaxLo","vmaxHi"))[,5]
mean_PFAB <- ft(df_test$beta, bcrit, E0, 
                para["ecrit"], V0, para["vmaxLo"], para["vmaxHi"], para["a"], para["b"])
var_PFAB <- dfdt(df_test$beta, bcrit, V0, para["vmaxLo"], para["vmaxHi"], para["a"], para["b"])

## Calculate MAE
MAE1 <- mean(abs(df_test$mean - syn_like_GP[[1]]))
MAE2 <- mean(abs(df_test$mean - syn_like_GP2[[1]]))
MAE3 <- mean(abs(df_test$mean - syn_like_GP3[[1]]))
MAE4 <- mean(abs(df_test$mean - mean_PFAB))
  
## Calculate RMSPE
RMSE1 <- sqrt(mean((df_test$mean - syn_like_GP[[1]])^2))
RMSE2 <- sqrt(mean((df_test$mean - syn_like_GP2[[1]])^2))
RMSE3 <- sqrt(mean((df_test$mean - syn_like_GP3[[1]])^2))
RMSE4 <- sqrt(mean((df_test$mean - mean_PFAB)^2))

## Calculate MAE
MAE1_v <- mean(abs(df_test$var - syn_like_GP[[2]]))
MAE2_v <- mean(abs(df_test$var - syn_like_GP2[[2]]))
MAE3_v <- mean(abs(df_test$var - syn_like_GP3[[2]]))
MAE4_v <- mean(abs(df_test$var - var_PFAB))

## Calculate RMSPE
RMSE1_v <- sqrt(mean((df_test$var - syn_like_GP[[2]])^2))
RMSE2_v <- sqrt(mean((df_test$var - syn_like_GP2[[2]])^2))
RMSE3_v <- sqrt(mean((df_test$var - syn_like_GP3[[2]])^2))
RMSE4_v <- sqrt(mean((df_test$var - var_PFAB)^2))

pred <- c(MAE1, MAE2, MAE3, MAE4, RMSE1, RMSE2, RMSE3, RMSE4,
          MAE1_v, MAE2_v, MAE3_v, MAE4_v, RMSE1_v, RMSE2_v, RMSE3_v, RMSE4_v)
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
