#' @export
synthetic_likefnc_GP_1D_stat <- function(object1, newdata, summary_stats, mean1){
  
  d1 <- object1
  mmat1 <- model.matrix(update(d1$f, NULL ~ .), newdata)
  mmat1_2 <- mmat1 + 1e-5
  
  ### mean, var for first summary stats
  
  s_in <- tf$constant(mmat1, dtype = "float64", name = "s")
  s_in_2 <- tf$constant(mmat1_2, dtype = "float64", name = "s")
  ndata <- nrow(d1$data)
  
  z_tf <- d1$z_tf_1

  obs_swarped <- d1$s_tf
  
  newdata_swarped <- s_in
  newdata_swarped_2 <- s_in_2

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
  
  pred_mean_1 <- d1$run(pred_mean_1) + mean1
  pred_mean_1_2 <- d1$run(pred_mean_1_2) + mean1
  
  pred_var_1 <- pmax((pred_mean_1_2 - pred_mean_1) / 1e-5, rep(1e-5, length(pred_mean_1)))
  pred_sd_1 <- sqrt(pred_var_1)
  
  z1 <- summary_stats

  syn_like <- - log(pred_sd_1) - 1/2 * (z1 - pred_mean_1)^2 / pred_var_1

  pred <- list(pred_mean_1, pred_var_1, syn_like)
  
}