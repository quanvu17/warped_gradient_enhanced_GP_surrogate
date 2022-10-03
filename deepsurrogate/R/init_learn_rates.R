#' @export
init_learn_rates <- function(sigma2y = 0.0005, covfun = 0.01, sigma2eta = 0.0001,
                             eta_mean = 0.1, eta_sd = 0.1, LFTpars = 0.01,
                             AFFpars = 0.01,
                             rho = 0.1) {
  
  list(sigma2y = sigma2y,
       covfun = covfun,
       sigma2eta = sigma2eta,
       eta_mean = eta_mean,
       eta_sd = eta_sd,
       LFTpars = LFTpars,
       AFFpars = AFFpars,
       rho = rho)
}
