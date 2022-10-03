tf_status <- new.env()

## Load TensorFlow and add the cholesky functions
.onLoad <- function(libname, pkgname) {
  
  tf <- reticulate::import("tensorflow", delay_load = TRUE)
  tf$cholesky_lower <- tf$linalg$cholesky
  tf$cholesky_upper <- function(x) tf$linalg$transpose(tf$linalg$cholesky(x))

}
