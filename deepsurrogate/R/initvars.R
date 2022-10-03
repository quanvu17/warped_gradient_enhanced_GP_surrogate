#' @export
initvars <- function(sigma2y = 0.1,
                     l_top_layer = 0.5,
                     sigma2eta_top_layer = 1,
                     nu = 1.5,
                     transeta_mean_init = list(AWU = -3,
                                               RBF = -0.8,
                                               LFT = 1,
                                               AFF_1D = 1,
                                               AFF_2D = 1),
                     transeta_mean_prior = list(AWU = -3,
                                                RBF = -0.8,
                                                LFT = NA),
                     transeta_sd_init = list(AWU = 0.01,
                                             RBF = 0.01,
                                             LFT = 0.01),
                     transeta_sd_prior = list(AWU = 2,
                                              RBF = 2,
                                              LFT = NA)) {

  list(sigma2y = sigma2y,
       sigma2eta_top_layer = sigma2eta_top_layer,
       l_top_layer = l_top_layer,
       nu = nu,
       transeta_mean_init = transeta_mean_init,
       transeta_mean_prior = transeta_mean_prior,
       transeta_sd_init =transeta_sd_init,
       transeta_sd_prior = transeta_sd_prior)
}
