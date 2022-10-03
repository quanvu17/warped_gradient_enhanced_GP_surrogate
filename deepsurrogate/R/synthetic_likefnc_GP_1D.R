#' @export
synthetic_likefnc_GP_1D <- function(object1, newdata, summary_stats, mean1){
  
  d1 <- object1
  mmat1 <- model.matrix(update(d1$f, NULL ~ .), newdata)
  mmat1_2 <- mmat1 + 1e-5
  
  ### mean, var for first summary stats
  
  s_in <- tf$constant(mmat1, dtype = "float64", name = "s")
  s_in_2 <- tf$constant(mmat1_2, dtype = "float64", name = "s")
  ndata <- nrow(d1$data)
  
  z_tf <- d1$z_tf_1
  s_in <- scale_0_5_tf(s_in, d1$scalings[[1]]$min, d1$scalings[[1]]$max)
  s_in_2 <- scale_0_5_tf(s_in_2, d1$scalings[[1]]$min, d1$scalings[[1]]$max)
  
  if (d1$method == "ML"){
  h_tf <- list(s_in)
  for (i in 1:d1$nlayers){
    h_tf[[i + 1]] <- d1$layers[[i]]$f(h_tf[[i]], d1$eta_tf[[i]]) %>%
        scale_0_5_tf(smin_tf = d1$scalings[[i + 1]]$min,
                        smax_tf = d1$scalings[[i + 1]]$max)
    } 
  
  h_tf_2 <- list(s_in_2)
  for (i in 1:d1$nlayers){
    h_tf_2[[i + 1]] <- d1$layers[[i]]$f(h_tf_2[[i]], d1$eta_tf[[i]]) %>%
        scale_0_5_tf(smin_tf = d1$scalings[[i + 1]]$min,
                        smax_tf = d1$scalings[[i + 1]]$max)
    } 
  
  obs_swarped <- d1$swarped_tf
  
  newdata_swarped <- h_tf[[d1$nlayers + 1]]
  newdata_swarped_2 <- h_tf_2[[d1$nlayers + 1]]

  K_obs <- cov_matern_tf_1D(x1 = obs_swarped, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  K_obs_star <- cov_matern_tf_1D(x1 = obs_swarped, x2 = newdata_swarped, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  K_obs_star_2 <- cov_matern_tf_1D(x1 = obs_swarped, x2 = newdata_swarped_2, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  
  K_obs_2 <- K_obs + d1$Sobs_tf_1
  
  Kobs_chol <- tf$cholesky(K_obs_2)
  Kobs_chol_z <- tf$matrix_solve(Kobs_chol, z_tf)
  Kobs_chol_z <- tf$matrix_solve(Kobs_chol, z_tf)
  Kobs_chol_star <- tf$matrix_solve(Kobs_chol, K_obs_star)
  Kobs_chol_star_2 <- tf$matrix_solve(Kobs_chol, K_obs_star_2)
  
  pred_mean_1 <- tf$matmul(tf$matrix_transpose(Kobs_chol_star), Kobs_chol_z)
  pred_mean_1_2 <- tf$matmul(tf$matrix_transpose(Kobs_chol_star_2), Kobs_chol_z)
  
  }
  
  pred_mean_1 <- d1$run(pred_mean_1) + mean1
  pred_mean_1_2 <- d1$run(pred_mean_1_2) + mean1
  
  pred_var_1 <- pmax((pred_mean_1_2 - pred_mean_1) / 1e-5, rep(1e-5, length(pred_mean_1)))
  pred_sd_1 <- sqrt(pred_var_1)
  
  z1 <- summary_stats

  syn_like <- - log(pred_sd_1) - 1/2 * (z1 - pred_mean_1)^2 / pred_var_1

  pred <- list(pred_mean_1, pred_var_1, syn_like)
  
}