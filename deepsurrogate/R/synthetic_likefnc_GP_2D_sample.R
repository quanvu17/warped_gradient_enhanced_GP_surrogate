synthetic_likefnc_GP_2D_sample <- function(object1, object2, newdata, summary_stats, mean1, mean2){
  
  d1 <- object1
  mmat1 <- model.matrix(update(d1$f, NULL ~ .), newdata)
  mmat1_2 <- mmat1 + cbind(rep(1e-5, nrow(mmat1)), rep(0, nrow(mmat1)))
  d2 <- object2
  mmat2 <- model.matrix(update(d2$f, NULL ~ .), newdata)
  mmat2_2 <- mmat2 + cbind(rep(1e-5, nrow(mmat2)), rep(0, nrow(mmat2)))
  
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
  
  K_obs <- cov_matern_tf_2D(x1 = obs_swarped, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  K_obs_star <- cov_matern_tf_2D(x1 = obs_swarped, x2 = newdata_swarped_0, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  K_star <- cov_matern_tf_2D(x1 = newdata_swarped_0, sigma2f = d1$sigma2_tf, alpha = tf$reciprocal(d1$l_tf))
  
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
  
  pred_mean_01 <- matrix(rep(0, 30 * nrow(newdata_swarped)), nrow = nrow(newdata_swarped))
  pred_var_01 <- matrix(rep(0, 30 * nrow(newdata_swarped)), nrow = nrow(newdata_swarped))

  for (i in 1:30){
    set.seed(i)
    normvec <- rnorm(length(pred_mean_1), 0, 1)
    pred_mean_i <- pred_GP_var_chol %*% normvec + pred_mean_1
    pred_mean_01[,i] <- pred_mean_i[1:nrow(newdata_swarped)]
    pred_mean_001 <- pred_mean_i[(nrow(newdata_swarped) + 1) : (2 * nrow(newdata_swarped))]
    pred_var_01[,i] <- (pred_mean_001 - pred_mean_01[,i]) / 1e-5
  }
  
  ### mean, var for second summary stats
  
  s_in <- tf$constant(mmat2, dtype = "float64", name = "s")
  s_in_2 <- tf$constant(mmat2_2, dtype = "float64", name = "s")
  ndata <- nrow(d2$data)
  
  z_tf <- d2$z_tf_1
  s_in <- scale_0_5_tf(s_in, d2$scalings[[1]]$min, d2$scalings[[1]]$max)
  s_in_2 <- scale_0_5_tf(s_in_2, d2$scalings[[1]]$min, d2$scalings[[1]]$max)
  
  if (d2$method == "ML"){
    h_tf <- list(s_in)
    for (i in 1:d2$nlayers){
      h_tf[[i + 1]] <- d2$layers[[i]]$f(h_tf[[i]], d2$eta_tf[[i]]) %>%
        scale_0_5_tf(smin_tf = d2$scalings[[i + 1]]$min,
                     smax_tf = d2$scalings[[i + 1]]$max)
    } 
    
    h_tf_2 <- list(s_in_2)
    for (i in 1:d2$nlayers){
      h_tf_2[[i + 1]] <- d2$layers[[i]]$f(h_tf_2[[i]], d2$eta_tf[[i]]) %>%
        scale_0_5_tf(smin_tf = d2$scalings[[i + 1]]$min,
                     smax_tf = d2$scalings[[i + 1]]$max)
    } 
    
    obs_swarped <- d2$swarped_tf
    
    newdata_swarped <- h_tf[[d2$nlayers + 1]]
    newdata_swarped_2 <- h_tf_2[[d2$nlayers + 1]]
    
    newdata_swarped_0 <- tf$concat(list(newdata_swarped, newdata_swarped_2), axis = 0L)
    
    K_obs <- cov_matern_tf_2D(x1 = obs_swarped, sigma2f = d2$sigma2_tf, alpha = tf$reciprocal(d2$l_tf))
    K_obs_star <- cov_matern_tf_2D(x1 = obs_swarped, x2 = newdata_swarped_0, sigma2f = d2$sigma2_tf, alpha = tf$reciprocal(d2$l_tf))
    K_star <- cov_matern_tf_2D(x1 = newdata_swarped_0, sigma2f = d2$sigma2_tf, alpha = tf$reciprocal(d2$l_tf))
    
    K_obs_2 <- K_obs + d2$Sobs_tf_1
    
    Kobs_chol <- tf$cholesky(K_obs_2)
    Kobs_chol_z <- tf$matrix_solve(Kobs_chol, z_tf)
    Kobs_chol_star <- tf$matrix_solve(Kobs_chol, K_obs_star)
    
    pred_mean_2 <- tf$matmul(tf$matrix_transpose(Kobs_chol_star), Kobs_chol_z)
    
    pred_GP_var <- K_star - tf$matmul(tf$matrix_transpose(Kobs_chol_star), Kobs_chol_star)
    pred_GP_var_chol <- tf$cholesky(pred_GP_var)
    
  }
  
  pred_mean_2 <- d2$run(pred_mean_2) + mean2
  pred_GP_var_chol <- d2$run(pred_GP_var_chol)
  
  pred_mean_02 <- matrix(rep(0, 30 * nrow(newdata_swarped)), nrow = nrow(newdata_swarped))
  pred_var_02 <- matrix(rep(0, 30 * nrow(newdata_swarped)), nrow = nrow(newdata_swarped))

  for (i in 1:30){
    set.seed(i)
    normvec <- rnorm(length(pred_mean_2), 0, 1)
    pred_mean_i <- pred_GP_var_chol %*% normvec + pred_mean_2
    pred_mean_02[,i] <- pred_mean_i[1:nrow(newdata_swarped)]
    pred_mean_002 <- pred_mean_i[(nrow(newdata_swarped) + 1) : (2 * nrow(newdata_swarped))]
    pred_var_02[,i] <- (pred_mean_002 - pred_mean_02[,i]) / 1e-5
  }
  
  z1 <- summary_stats[1]
  z2 <- summary_stats[2]
  
  syn_like_0 <- matrix(rep(0, 30 * nrow(newdata_swarped)), nrow = nrow(newdata_swarped))
  for (i in 1:30){
    Cost1 <- - 1/2 * log(pred_var_01[,i]) - 1/2 * (z1 - pred_mean_01[,i])^2 / pred_var_01[,i]
    Cost2 <- - 1/2 * log(pred_var_02[,i]) - 1/2 * (z2 - pred_mean_02[,i])^2 / pred_var_02[,i]
    syn_like_0[,i] <- Cost1 + Cost2
  }

  pred <- list(pred_mean_01, pred_mean_02, pred_var_01, pred_var_02, syn_like_0)
  
}