## Matern Covariance
cov_matern_tf_1D <- function(x1, x2 = x1, sigma2f, alpha) {

  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
    x11 <- x1[, 1, drop = FALSE]
    x21 <- x2[, 1, drop = FALSE]
    D <- x11 - tf$matrix_transpose(x21)

    D_abs <- tf$abs(D)
    D_abs <- D_abs + 1e-10
  
  aD <- alpha * D_abs

  K <- tf$multiply(sigma2f, tf$multiply((tf$constant(1, dtype = "float64") + aD), tf$exp(-aD)))
  return(K)
}

d_cov_matern_tf_1D <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D <- x11 - tf$matrix_transpose(x21)
  
  D_abs <- tf$abs(D)
  D_abs <- D_abs + 1e-10
  
  aD <- alpha * D_abs
  
  dK <- tf$multiply(-sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), D))
  return(dK)
}

d2_cov_matern_tf_1D <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D <- x11 - tf$matrix_transpose(x21)
  
  D_abs <- tf$abs(D)
  D_abs <- D_abs + 1e-10
  
  aD <- alpha * D_abs
  
  d2K <- tf$multiply(sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), tf$constant(1, dtype = "float64") - alpha * tf$square(D) / D_abs))
  return(d2K)
}

cov_matern_tf_2D <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D1 <- x11 - tf$matrix_transpose(x21)
  
  x12 <- x1[, 2, drop = FALSE]
  x22 <- x2[, 2, drop = FALSE]
  D2 <- x12 - tf$matrix_transpose(x22)
  
  D <- tf$square(D1) + tf$square(D2)
  
  D_abs <- tf$sqrt(D + 1e-10)
  
  aD <- alpha * D_abs
  
  K <- tf$multiply(sigma2f, tf$multiply((tf$constant(1, dtype = "float64") + aD), tf$exp(-aD)))
  return(K)
}

d_cov_matern_tf_2D_v1 <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D1 <- x11 - tf$matrix_transpose(x21)
  
  x12 <- x1[, 2, drop = FALSE]
  x22 <- x2[, 2, drop = FALSE]
  D2 <- x12 - tf$matrix_transpose(x22)
  
  D <- tf$square(D1) + tf$square(D2)
  D_abs <- tf$sqrt(D + 1e-10)
  
  aD <- alpha * D_abs
  
  dK <- tf$multiply(-sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), D1))
  return(dK)
}

d_cov_matern_tf_2D_v2 <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D1 <- x11 - tf$matrix_transpose(x21)
  
  x12 <- x1[, 2, drop = FALSE]
  x22 <- x2[, 2, drop = FALSE]
  D2 <- x12 - tf$matrix_transpose(x22)
  
  D <- tf$square(D1) + tf$square(D2)
  D_abs <- tf$sqrt(D + 1e-10)
  
  aD <- alpha * D_abs
  
  dK <- tf$multiply(-sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), D2))
  return(dK)
}

d2_cov_matern_tf_2D_v1 <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D1 <- x11 - tf$matrix_transpose(x21)
  
  x12 <- x1[, 2, drop = FALSE]
  x22 <- x2[, 2, drop = FALSE]
  D2 <- x12 - tf$matrix_transpose(x22)
  
  D <- tf$square(D1) + tf$square(D2)
  D_abs <- tf$sqrt(D + 1e-10)
  
  aD <- alpha * D_abs
  
  d2K <- tf$multiply(sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), tf$constant(1, dtype = "float64") - alpha * tf$square(D1) / D_abs))
  return(d2K)
}

d2_cov_matern_tf_2D_v2 <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D1 <- x11 - tf$matrix_transpose(x21)
  
  x12 <- x1[, 2, drop = FALSE]
  x22 <- x2[, 2, drop = FALSE]
  D2 <- x12 - tf$matrix_transpose(x22)
  
  D <- tf$square(D1) + tf$square(D2)
  D_abs <- tf$sqrt(D + 1e-10)
  
  aD <- alpha * D_abs
  
  d2K <- tf$multiply(sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), tf$constant(1, dtype = "float64") - alpha * tf$square(D2) / D_abs))
  return(d2K)
}

d2_cov_matern_tf_2D_v12 <- function(x1, x2 = x1, sigma2f, alpha) {
  
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  D <- tf$constant(matrix(0, n1, n2), name='D', dtype = tf$float64)
  
  x11 <- x1[, 1, drop = FALSE]
  x21 <- x2[, 1, drop = FALSE]
  D1 <- x11 - tf$matrix_transpose(x21)
  
  x12 <- x1[, 2, drop = FALSE]
  x22 <- x2[, 2, drop = FALSE]
  D2 <- x12 - tf$matrix_transpose(x22)
  
  D <- tf$square(D1) + tf$square(D2)
  D_abs <- tf$sqrt(D + 1e-10)
  
  aD <- alpha * D_abs
  
  d2K <- (-tf$constant(1, dtype = "float64")) * tf$multiply(sigma2f * tf$square(alpha), tf$multiply(tf$exp(-aD), D1 * D2 / D_abs))
  return(d2K)
}