#' @export
surrogate_GP <- function(f, data, method = c("ML", "VB"), layers = layers,
                           model = c("meanonly", "meanvar", "meansd"),
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
  Sobs_tf_2 <- tf$diag(tf$constant(noise_var_v + 1e-3, dtype = "float64"))

  Sobs_tf <- tf$concat(list(tf$concat(list(Sobs_tf_1, tf$zeros(shape=c(ndata,ndata), dtype=tf$float64)), axis=1L),
                            tf$concat(list(tf$zeros(shape=c(ndata,ndata), dtype=tf$float64), Sobs_tf_2), axis=1L)),axis=0L)
  
  ## Prior variance of the process
  sigma2 <- var(data[[depvar1]]) #par_init$sigma2eta_top_layer
  logsigma2_tf <- tf$Variable(log(sigma2), name = "sigma2eta", dtype = "float64")
  sigma2_tf <- tf$exp(logsigma2_tf)
  
  ## Length scale of process
  l <- par_init$l_top_layer
  logl_tf <-tf$Variable(matrix(log(l)), name = "l", dtype = "float64")
  l_tf <- tf$exp(logl_tf)
  
  nlayers <- length(layers)
  scalings <- list(scale_lims_tf(s_tf))
  s_tf <- scale_0_5_tf(s_tf, scalings[[1]]$min, scalings[[1]]$max)
  
  if(method == "ML") {
    ## Do the warping
    transeta_tf <- eta_tf <- swarped_tf <- list()
    swarped_tf[[1]] <- s_tf
    layer_type <- "AWU"
    
    swarped_final_d <- tf$reshape(tf$constant(rep(1, ndata), dtype = "float64"), c(ndata, 1L))
    swarped_final_d2 <- tf$reshape(tf$constant(rep(0, ndata), dtype = "float64"), c(ndata, 1L))
    
    for(i in 1:nlayers) {
      transeta_tf[[i]] <- tf$Variable(matrix(rep(par_init$transeta_mean_init[[layer_type]], layers[[i]]$r)),
                                      name = paste0("eta", i), dtype = "float64")
      
      eta_tf[[i]] <- layers[[i]]$trans(transeta_tf[[i]]) # ensure positivity for some variables
      
      swarped_tf[[i + 1]] <- layers[[i]]$f(swarped_tf[[i]], eta_tf[[i]])
      
        scalings[[i + 1]] <- scale_lims_tf(swarped_tf[[i + 1]])
        swarped_tf[[i + 1]] <- scale_0_5_tf(swarped_tf[[i + 1]], scalings[[i + 1]]$min, scalings[[i + 1]]$max)
      
      swarped_new_d <- layers[[i]]$f_d(swarped_tf[[i]], eta_tf[[i]]) * 1/(scalings[[i]]$max[,1] - scalings[[i]]$min[,1])
      swarped_new_d2 <- layers[[i]]$f_d2(swarped_tf[[i]], eta_tf[[i]]) * 1/(scalings[[i]]$max[,1] - scalings[[i]]$min[,1])
      
      swarped_old_d <- swarped_final_d
      swarped_old_d2 <- swarped_final_d2
      
      swarped_final_d <- swarped_old_d * swarped_new_d
      swarped_final_d2 <- swarped_old_d2 * swarped_new_d + swarped_new_d2 * (swarped_old_d)^2
    }
    
    swarped_final <- swarped_tf[[nlayers + 1]]

  } 
  
  ##############################################################
  ##Training
  if (model == "meanonly"){
    if(method == "ML") {
      NMLL <- logmarglik_1D_warp(s_in = swarped_final,
                                Sobs_tf = Sobs_tf_1,
                                l_tf = l_tf,
                                sigma2_tf = sigma2_tf,
                                z_tf = z_tf_1)
      Cost <- NMLL$Cost
    }
  }
  if (model == "meanvar"){
  if(method == "ML") {
    NMLL <- logmarglik_d_1D_warp(s_in = swarped_final,
                                 Sobs_tf = Sobs_tf,
                                 l_tf = l_tf,
                                 sigma2_tf = sigma2_tf,
                                 z_tf = z_tf,
                                 swarped_d = swarped_final_d,
                                 swarped_d2 = swarped_final_d2)
    Cost <- NMLL$Cost
  }
  }
  
  ## Optimisers for top layer
  # trains2y = (tf$train$GradientDescentOptimizer(learn_rates$sigma2y))$minimize(Cost, var_list = list(logsigma2y_tf_1, logsigma2y_tf_2))
  traincovfun = (tf$train$AdamOptimizer(learn_rates$covfun))$minimize(Cost, var_list = list(logl_tf))
  trains2eta = (tf$train$AdamOptimizer(learn_rates$sigma2eta))$minimize(Cost, var_list = list(logsigma2_tf))
  
  ## Optimisers for eta (all hidden layers except LFT)
  nLFTlayers <- sum(sapply(layers, function(l) l$name) == "LFT")
  LFTidx <- which(sapply(layers, function(l) l$name) == "LFT")
  notLFTidx <- setdiff(1:nlayers, LFTidx)
  opt_eta <- (nlayers > 0) & (nLFTlayers < nlayers)
  if(opt_eta)
    if(method == "ML") {
      traineta_mean <- (tf$train$AdamOptimizer(learn_rates$eta_mean))$minimize(Cost, var_list = list(transeta_tf[notLFTidx]))
      #traineta_mean_1 <- (tf$train$AdamOptimizer(1e-3))$minimize(Cost, var_list = list(logl_tf_0))
    }
  
  init <- tf$global_variables_initializer()
  run <- tf$Session()$run
  run(init)
  Objective <- rep(0, 2*nsteps)
  
  if(method == "ML") {
    negcostname <- "Likelihood"
  } 
  
  cat("Learning weight parameters... \n")
  for(i in 1:nsteps) {
    if(opt_eta) run(traineta_mean)
    if(nLFTlayers > 0) run(trainLFTpars)
    thisML <- -run(Cost)
    if((i %% 10) == 0)
      cat(paste0("Step ", i, " ... ", negcostname, ": ", thisML, "\n"))
    Objective[i] <- thisML
  }
  
 #cat("Measurement-error variance and cov. fn. parameters... \n")
 #for(i in (nsteps + 1):(2 * nsteps)) {
 #  if(nLFTlayers > 0) run(trainLFTpars)
 #  # run(trains2y)
 #  run(traincovfun)
 #  run(trains2eta)
 #  thisML <- -run(Cost)
 #  if((i %% 10) == 0)
 #    cat(paste("Step ", i, " ... ", negcostname, ": ", thisML, "\n"))
 #  Objective[i] <- thisML
 #}
  
  cat("Updating everything... \n")
  for(i in (nsteps + 1):(2 * nsteps)) {
    if(opt_eta) run(traineta_mean)
    if(nLFTlayers > 0) run(trainLFTpars)
    # run(trains2y)
    run(traincovfun)
    run(trains2eta)
    thisML <- -run(Cost)
    if((i %% 10) == 0)
      cat(paste0("Step ", i, " ... ", negcostname, ": ", thisML, "\n"))
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
                       swarped_tf = swarped_tf[[nlayers + 1]],
                       layers = layers,
                       nlayers = nlayers,
                       eta_tf = eta_tf,
                       sigma2_tf = sigma2_tf,
                       l_tf = l_tf,
                       Sobs_tf_1 = Sobs_tf_1,
                       scalings = scalings
  )
  
  class(deepspat.obj) <- "surrogate_GP"
  deepspat.obj
}
