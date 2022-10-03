#' @export
AWU <- function(r = 50L, dim = 1L, grad = 200, lims = c(-0.5, 0.5)) {

  ## Parameters appearing in sigmoid (grad, loc)
  theta <- matrix(c(grad, 0), nrow = r - 1, ncol = 2, byrow = TRUE)
  theta[, 2] <- seq(lims[1], lims[2], length.out = (r - 1) + 2)[-c(1, (r - 1) + 2)]

  theta_steep_unclipped_tf <- tf$constant(theta[, 1, drop = FALSE], name = "thetasteep", dtype = "float64")
  theta_steep_tf <- tf$clip_by_value(theta_steep_unclipped_tf, 0, 200)
  theta_locs_tf <- tf$constant(theta[, 2, drop = FALSE], name = "thetalocs", dtype = "float64")
  theta_tf <- tf$concat(list(theta_steep_tf, theta_locs_tf), 1L)

  f = function(s_tf, eta_tf) {
    PHI_tf <- tf$concat(list(s_tf[, dim, drop = FALSE],
                  sigmoid_tf(s_tf[, dim, drop = FALSE], theta_tf)), 1L)
    swarped <-  tf$matmul(PHI_tf, eta_tf)
    slist <- lapply(1:ncol(s_tf), function(i) s_tf[, i, drop = FALSE])
    slist[[dim]] <- swarped
    sout_tf <- tf$concat(slist, axis = 1L)
  }
  
  
  f_d = function(s_tf, eta_tf) {
    if (dim == 1L){
    PHI_tf <- tf$concat(list(tf$constant(1, shape = c(nrow(s_tf), 1L), dtype = "float64"),
                             sigmoid_d_tf(s_tf[, dim, drop = FALSE], theta_tf)), 1L)
    swarped_d <-  tf$matmul(PHI_tf, eta_tf)
    }
    else{
    swarped_d <- tf$constant(1, shape = c(nrow(s_tf), 1L), dtype = "float64")
    }
    swarped_d
  }
  
  
  f_d2 = function(s_tf, eta_tf) {
    if (dim == 1L){
    PHI_tf <- tf$concat(list(tf$constant(0, shape = c(nrow(s_tf), 1L), dtype = "float64"),
                             sigmoid_d2_tf(s_tf[, dim, drop = FALSE], theta_tf)), 1L)
    swarped_d2 <-  tf$matmul(PHI_tf, eta_tf)
    }
    else{
    swarped_d2 <- tf$constant(0, shape = c(nrow(s_tf), 1L), dtype = "float64")
    }
    swarped_d2
  }
  

  fMC = function(s_tf, eta_tf) {
    PHI_tf <- list(s_tf[, , dim, drop = FALSE],
                   sigmoid_tf(s_tf[, , dim, drop = FALSE], theta_tf)) %>%
      tf$concat(2L)
    swarped <-  tf$matmul(PHI_tf, eta_tf)
    slist <- lapply(1:ncol(s_tf[1, , ]), function(i) s_tf[, , i, drop = FALSE])
    slist[[dim]] <- swarped
    sout_tf <- tf$concat(slist, axis = 2L)
  }

  f_dMC = function(s_tf, eta_tf) {
    if (dim == 1L){
    PHI_tf <- tf$concat(list(tf$constant(1, shape = c(nrow(s_tf), ncol(s_tf), 1L), dtype = "float64"),
                             sigmoid_dMC_tf(s_tf[, , dim, drop = FALSE], theta_tf)), 2L)
    swarped_d <-  tf$matmul(PHI_tf, eta_tf)
    }
    else{
    swarped_d  <- tf$constant(1, shape = c(nrow(s_tf), ncol(s_tf), 1L), dtype = "float64")
    }
  }
  
  fR = function(s, eta) {
    PHI_list <- list(s[, dim, drop = FALSE],
                     sigmoid(s[, dim, drop = FALSE], theta))
    PHI <- do.call("cbind", PHI_list)
    swarped <-  PHI %*% eta
    slist <- lapply(1:ncol(s), function(i) s[, i, drop = FALSE])
    slist[[dim]] <- swarped
    sout <- do.call("cbind", slist)
  }

  list(list(f = f,
            f_d = f_d,
            f_d2 = f_d2,
            fR = fR,
            fMC = fMC,
            f_dMC = f_dMC,
            r = r,
            trans = tf$exp,
            fix_weights = FALSE,
            name = "AWU"))
}
