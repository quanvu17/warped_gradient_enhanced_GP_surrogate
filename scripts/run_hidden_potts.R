cat("Fitting the surrogate model... \n")
source("scripts/potts_02_surrogate_fit.R")

cat("Running delayed-acceptance MCMC... \n")
source("scripts/hidden_potts_01_inference_delayaccept.R")

cat("Running exchange algorithm... \n")
source("scripts/hidden_potts_01_inference_exchange.R")

cat("Running surrogate MCMC... \n")
source("scripts/hidden_potts_01_inference_surrogate.R")

cat("Summarizing results... \n")
source("scripts/hidden_potts_02_summary.R")
