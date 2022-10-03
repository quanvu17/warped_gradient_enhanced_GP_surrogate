#' @export
surrogate_GP_stat <- function(f, data, method = c("ML", "VB"),
                           par_init = initvars(),
                           learn_rates = init_learn_rates(),
                           MC = 10L, nsteps,
                           noise_var_m, noise_var_v) {
  
  stopifnot(is(f, "formula"))
  stopifnot(is(data, "data.frame"))
  method = match.arg(method)
  mmat <- model.matrix(f, data)
  
  s_tf <- tf$constant(mmat, name = "s", dtype = "float64")
  
  depvar <- get_depvars_multivar(f)
  depvar1 <- depvar[1]
  depvar2 <- depvar[2]
  z_tf_1 <- tf$constant(as.matrix(data[[depvar1]]), name = 'z1', dtype = 'float64')
  z_tf_2 <- tf$constant(as.matrix(data[[depvar2]]), name = 'z2', dtype = 'float64')
  z_tf <- tf$concat(list(z_tf_1, z_tf_2), axis=0L)
  ndata <- nrow(data)

  ## Error
  Sobs_tf_1 <- tf$diag(tf$constant(noise_var_m + 1e-3, dtype = "float64"))

  ## Prior variance of the process
  sigma2 <- var(data[[depvar1]]) #par_init$sigma2eta_top_layer
  logsigma2_tf <- tf$Variable(log(sigma2), name = "sigma2eta", dtype = "float64")
  sigma2_tf <- tf$exp(logsigma2_tf)
  
  ## Length scale of process
  l <- par_init$l_top_layer
  logl_tf <-tf$Variable(matrix(log(l)), name = "l", dtype = "float64")
  l_tf <- tf$exp(logl_tf)
  
  ##############################################################
  ##Training
    if(method == "ML") {
      NMLL <- logmarglik_1D_warp(s_in = s_tf,
                                Sobs_tf = Sobs_tf_1,
                                l_tf = l_tf,
                                sigma2_tf = sigma2_tf,
                                z_tf = z_tf_1)
      Cost <- NMLL$Cost
    }
 
  ## Optimisers for top layer
  # trains2y = (tf$train$GradientDescentOptimizer(learn_rates$sigma2y))$minimize(Cost, var_list = list(logsigma2y_tf_1, logsigma2y_tf_2))
  traincovfun = (tf$train$AdamOptimizer(learn_rates$covfun))$minimize(Cost, var_list = list(logl_tf))
  trains2eta = (tf$train$AdamOptimizer(learn_rates$sigma2eta))$minimize(Cost, var_list = list(logsigma2_tf))
  
  init <- tf$global_variables_initializer()
  run <- tf$Session()$run
  run(init)
  Objective <- rep(0, 2*nsteps)
  
  if(method == "ML") {
    negcostname <- "Likelihood"
  } 
  
  cat("Measurement-error variance and cov. fn. parameters... \n")
  for(i in 1:(2 * nsteps)) {
    # run(trains2y)
    run(traincovfun)
    run(trains2eta)
    thisML <- -run(Cost)
    if((i %% 10) == 0)
      cat(paste("Step ", i, " ... ", negcostname, ": ", thisML, "\n"))
    Objective[i] <- thisML
  }
  
  deepspat.obj <- list(f = f,
                       data = data,
                       Cost = Cost,
                       method = method,
                       MC = MC,
                       run = run,
                       f = f,
                       s_tf = s_tf,
                       negcost = Objective,
                       z_tf_1 = z_tf_1,
                       z_tf_2 = z_tf_2,
                       sigma2_tf = sigma2_tf,
                       l_tf = l_tf,
                       Sobs_tf_1 = Sobs_tf_1
  )
  
  class(deepspat.obj) <- "surrogate_GP_stat"
  deepspat.obj
}
