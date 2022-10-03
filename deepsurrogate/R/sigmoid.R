sigmoid <- function(x, theta) {
  PHI <- list()
  for(i in 1:nrow(theta)) {
    PHI[[i]] <- 1 / (1 + exp(-theta[i, 1] * (x - theta[i, 2])))
  }
  do.call("cbind", PHI)
}

sigmoid_tf <- function(x, theta) {

  theta1 <- tf$transpose(theta[, 1, drop = FALSE])
  theta2 <- tf$transpose(theta[, 2, drop = FALSE])

  tf$subtract(x, theta2) %>%
    tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
    tf$multiply(theta1) %>%
    tf$exp() %>%
    tf$add(tf$constant(1L, dtype = "float64")) %>%
    tf$reciprocal()
}

sigmoid_d_tf <- function(x, theta) {
  
  theta1 <- tf$transpose(theta[, 1, drop = FALSE])
  theta2 <- tf$transpose(theta[, 2, drop = FALSE])
  
  E <- tf$subtract(x, theta2) %>%
    tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
    tf$multiply(theta1) %>%
    tf$exp()
  
   tf$multiply(theta1, E) / tf$square(tf$constant(1, dtype = "float64") + E)
  
  #D <- tf$subtract(x, theta2) %>%
  #  tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
  #  tf$multiply(theta1) %>%
  #  tf$exp() %>%
  #  tf$add(tf$constant(1L, dtype = "float64")) %>%
  #  tf$square() %>%
  #  tf$reciprocal()
  #
  #tf$subtract(x, theta2) %>%
  #  tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
  #  tf$multiply(theta1) %>%
  #  tf$exp() %>%
  #  tf$multiply(theta1) %>%
  #  tf$multiply(D)
}

sigmoid_dMC_tf <- function(x, theta) {
  theta1 <- tf$transpose(theta[, 1, drop = FALSE])
  theta1 <- theta1 %>% tf$reshape(c(1L, 1L, ncol(theta1))) %>% tf$tile(c(nrow(x), 1L, 1L))
  theta2 <- tf$transpose(theta[, 2, drop = FALSE])
  theta2 <- theta2 %>% tf$reshape(c(1L, 1L, ncol(theta2))) %>% tf$tile(c(nrow(x), 1L, 1L))
  
  E <- tf$subtract(x, theta2) %>%
    tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
    tf$multiply(theta1) %>%
    tf$exp()
  
  tf$multiply(theta1, E) / tf$square(tf$constant(1, dtype = "float64") + E)
}

sigmoid_d2_tf <- function(x, theta) {
  
  theta1 <- tf$transpose(theta[, 1, drop = FALSE])
  theta2 <- tf$transpose(theta[, 2, drop = FALSE])
  
  E <- tf$subtract(x, theta2) %>%
    tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
    tf$multiply(theta1) %>%
    tf$exp()
  
   (tf$constant(2, dtype = "float64") * tf$square(theta1) * tf$square(E) - tf$square(theta1) * E * (tf$constant(1, dtype = "float64") + E) ) / tf$pow(tf$constant(1, dtype = "float64") + E, tf$constant(3, dtype = "float64"))
  
 #D <- tf$subtract(x, theta2) %>%
 #  tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
 #  tf$multiply(theta1) %>%
 #  tf$exp() %>%
 #  tf$add(tf$constant(1L, dtype = "float64")) %>%
 #  tf$pow(tf$constant(3, dtype = "float64")) %>%
 #  tf$reciprocal()
 #
 #A <- tf$subtract(x, theta2) %>%
 #  tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
 #  tf$multiply(theta1) %>%
 #  tf$exp() %>%
 #  tf$square() %>%
 #  tf$multiply(tf$square(theta1)) %>%
 #  tf$multiply(D)
 # 
 #B <- tf$subtract(x, theta2) %>%
 #  tf$multiply(tf$constant(-1L, dtype = "float64")) %>%
 #  tf$multiply(theta1) %>%
 #  tf$exp() %>%
 #  tf$multiply(tf$square(theta1)) %>%
 #  tf$multiply(D)
 # 
 #A - B
}