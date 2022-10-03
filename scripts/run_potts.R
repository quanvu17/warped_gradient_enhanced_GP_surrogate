cat("Simulating a data set from the true parameter... \n")
source("scripts/potts_00_simulate_data.R")

cat("Simulating the pseudo-data to fit the surrogate model... \n")
source("scripts/potts_01_sample.R")

cat("Fitting the surrogate model... \n")
source("scripts/potts_02_surrogate_fit.R")

cat("Comparing the different surrogate models... \n")
source("scripts/potts_03_validation.R")

cat("Running delayed-acceptance MCMC... \n")
source("scripts/potts_04_inference_delayaccept.R")

cat("Running exchange algorithm... \n")
source("scripts/potts_04_inference_exchange.R")

cat("Running importance sampling... \n")
source("scripts/potts_04_inference_importance.R")

cat("Summarizing results... \n")
source("scripts/potts_05_summary.R")
