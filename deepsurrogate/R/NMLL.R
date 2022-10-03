logmarglik_1D_warp <- function(s_in, Sobs_tf, l_tf, sigma2_tf, z_tf) {
  
  n <- nrow(s_in)
  
  C11_tf <-  cov_matern_tf_1D(x1 = s_in,
                              sigma2f = sigma2_tf,
                              alpha = tf$reciprocal(l_tf)
  )
  
  SY_tf <- C11_tf
  
  ## Compute posterior distribution of weights and the Cholesky factor
  SZ_tf <- SY_tf + Sobs_tf
  L_tf <- tf$cholesky_lower(SZ_tf)
  a <- tf$matrix_solve(L_tf, z_tf)
  
  Part1 <- tf$constant(-0.5, dtype = "float64") * tf$constant(2, dtype = "float64") * tf$reduce_sum(tf$log(tf$linalg$diag_part(L_tf)))
  Part2 <- tf$constant(-0.5, dtype = "float64") * tf$reduce_sum(tf$square(a))

  Cost <- -(Part1 + Part2)
  
  list(Cost = Cost)
  
}

logmarglik_d_1D_warp <- function(s_in, Sobs_tf, l_tf, sigma2_tf, z_tf, swarped_d, swarped_d2) {
  
  n <- nrow(s_in)
  
  C11_tf <-  cov_matern_tf_1D(x1 = s_in,
                                    sigma2f = sigma2_tf,
                                    alpha = tf$reciprocal(l_tf)
  )
  
  C21_tf0 <- d_cov_matern_tf_1D(x1 = s_in,
                           sigma2f = sigma2_tf,
                           alpha = tf$reciprocal(l_tf)
  )
  
  C21_tf <-  tf$multiply(swarped_d, C21_tf0)
  
  C12_tf <-  tf$transpose(C21_tf)
  
  C22_tf0 <-  d2_cov_matern_tf_1D(x1 = s_in,
                                       sigma2f = sigma2_tf,
                                       alpha = tf$reciprocal(l_tf)
  )
  
  C22_tf00 <- tf$multiply(swarped_d2, tf$reshape(tf$diag_part(C21_tf0), c(n, 1L)))[,1]
  
  C22_tf <- tf$multiply(tf$transpose(swarped_d), tf$multiply(swarped_d, C22_tf0)) + tf$diag(C22_tf00)
  
  SY_tf <- tf$concat(list(tf$concat(list(C11_tf,C12_tf), axis=1L),
                          tf$concat(list(C21_tf,C22_tf), axis=1L)), axis=0L)
  
  ## Compute posterior distribution of weights and the Cholesky factor
  SZ_tf <- SY_tf + Sobs_tf
  L_tf <- tf$cholesky_lower(SZ_tf)
  a <- tf$matrix_solve(L_tf, z_tf)
  
  Part1 <- tf$constant(-0.5, dtype = "float64") * tf$constant(2, dtype = "float64") * tf$reduce_sum(tf$log(tf$linalg$diag_part(L_tf)))
  Part2 <- tf$constant(-0.5, dtype = "float64") * tf$reduce_sum(tf$square(a))

  Cost <- -(Part1 + Part2)
  
  list(Cost = Cost)
  
}

logmarglik_2D_warp <- function(s_in, Sobs_tf, l_tf, sigma2_tf, z_tf) {
  
  n <- nrow(s_in)
  
  C11_tf <-  cov_matern_tf_2D(x1 = s_in,
                              sigma2f = sigma2_tf,
                              alpha = tf$reciprocal(l_tf)
  )
  
  SY_tf <- C11_tf
  
  ## Compute posterior distribution of weights and the Cholesky factor
  SZ_tf <- SY_tf + Sobs_tf
  L_tf <- tf$cholesky_lower(SZ_tf)
  a <- tf$matrix_solve(L_tf, z_tf)
  
  Part1 <- tf$constant(-0.5, dtype = "float64") * tf$constant(2, dtype = "float64") * tf$reduce_sum(tf$log(tf$linalg$diag_part(L_tf)))
  Part2 <- tf$constant(-0.5, dtype = "float64") * tf$reduce_sum(tf$square(a))
  
  Cost <- -(Part1 + Part2)
  
  list(Cost = Cost)
  
}

logmarglik_d_2D1_warp <- function(s_in, Sobs_tf, l_tf, sigma2_tf, z_tf,
                                  swarped_d_v1, swarped_d2_v1) {
  
  n <- nrow(s_in)
  
  C11_tf <-  cov_matern_tf_2D(x1 = s_in,
                              sigma2f = sigma2_tf,
                              alpha = tf$reciprocal(l_tf)
  )
  
  C21_tf0 <- d_cov_matern_tf_2D_v1(x1 = s_in,
                                   sigma2f = sigma2_tf,
                                   alpha = tf$reciprocal(l_tf)
  )
  
  C21_tf <-  tf$multiply(swarped_d_v1, C21_tf0)
  
  C12_tf <-  tf$transpose(C21_tf)
  
  C22_tf0 <-  d2_cov_matern_tf_2D_v1(x1 = s_in,
                                     sigma2f = sigma2_tf,
                                     alpha = tf$reciprocal(l_tf)
  )
  
  C22_tf00 <- tf$multiply(swarped_d2_v1, tf$reshape(tf$diag_part(C21_tf0), c(n, 1L)))[,1]
  
  C22_tf <- tf$multiply(tf$transpose(swarped_d_v1), tf$multiply(swarped_d_v1, C22_tf0)) + tf$diag(C22_tf00)
  
  SY_tf <- tf$concat(list(tf$concat(list(C11_tf,C12_tf), axis=1L),
                          tf$concat(list(C21_tf,C22_tf), axis=1L)), axis=0L)
  
  ## Compute posterior distribution of weights and the Cholesky factor
  SZ_tf <- SY_tf + Sobs_tf
  L_tf <- tf$cholesky_lower(SZ_tf)
  #a <- tf$matrix_inverse(SZ_tf)
  a <- tf$matrix_solve(L_tf, z_tf)
  
  Part1 <- tf$constant(-0.5, dtype = "float64") * tf$constant(2, dtype = "float64") * tf$reduce_sum(tf$log(tf$linalg$diag_part(L_tf)))
  Part2 <- tf$constant(-0.5, dtype = "float64") * tf$reduce_sum(tf$square(a))
  #Part2 <- -0.5 * tf$reduce_sum(tf$multiply(tf$multiply(tf$matrix_transpose(z_tf), a), z_tf))
  
  Cost <- -(Part1 + Part2)
  
  list(Cost = Cost)
  
}