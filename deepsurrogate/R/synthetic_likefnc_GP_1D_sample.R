#' @export
synthetic_likefnc_GP_1D_sample <- function(object1, newdata, summary_stats, mean1){
  
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

  newdata_swarped_0 <- tf$concat(list(newdata_swarped, newdata_swarped_2), axis = 0L)
  
  K_obs <- cov_matern_tf_1D(x1 = obs_swarped, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  K_obs_star <- cov_matern_tf_1D(x1 = obs_swarped, x2 = newdata_swarped_0, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  K_star <- cov_matern_tf_1D(x1 = newdata_swarped_0, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  
  K_obs_2 <- K_obs + d1$Sobs_tf_1
  
  Kobs_chol <- tf$cholesky(K_obs_2)
  Kobs_chol_z <- tf$matrix_solve(Kobs_chol, z_tf)
  Kobs_chol_star <- tf$matrix_solve(Kobs_chol, K_obs_star)

  pred_mean_1 <- tf$matmul(tf$matrix_transpose(Kobs_chol_star), Kobs_chol_z)

  pred_GP_var <- K_star - tf$matmul(tf$matrix_transpose(Kobs_chol_star), Kobs_chol_star)
  pred_GP_var_chol <- tf$cholesky(pred_GP_var)
  
  }
  
  pred_mean_1 <- d1$run(pred_mean_1) + mean1
  pred_GP_var_chol <- d1$run(pred_GP_var_chol)
  
  pred_mean_0 <- matrix(rep(0, 30 * length(newdata_swarped)), nrow = length(newdata_swarped))
  pred_var_0 <- matrix(rep(0, 30 * length(newdata_swarped)), nrow = length(newdata_swarped))
  syn_like_0 <- matrix(rep(0, 30 * length(newdata_swarped)), nrow = length(newdata_swarped))
  z1 <- summary_stats
  
  for (i in 1:30){
    set.seed(i)
    normvec <- rnorm(length(pred_mean_1), 0, 1)
    pred_mean_i <- pred_GP_var_chol %*% normvec + pred_mean_1
    pred_mean_0[,i] <- pred_mean_i[1:length(newdata_swarped)]
    pred_mean_00 <- pred_mean_i[(length(newdata_swarped) + 1) : (2 * length(newdata_swarped))]
    pred_var_0[,i] <- (pred_mean_00 - pred_mean_0[,i]) / 1e-5
    pred_sd_0 <- sqrt(pred_var_0[,i])
    syn_like_0[,i] <- - log(pred_sd_0) - 1/2 * (z1 - pred_mean_0[,i])^2 / pred_var_0[,i]
  }
  
  pred <- list(pred_mean_1, pred_mean_0, pred_var_0, syn_like_0)
  
}