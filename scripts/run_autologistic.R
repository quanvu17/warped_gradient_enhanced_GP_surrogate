cat("Simulating the pseudo-data to fit the surrogate model... \n")
source("scripts/autologistic_01_sample.R")

cat("Fitting the surrogate model... \n")
source("scripts/autologistic_02_surrogate_fit.R")

cat("Running exchange algorithm... \n")
source("scripts/autologistic_03_inference_exchange.R")

cat("Running importance sampling... \n")
source("scripts/autologistic_03_inference_importance.R")

cat("Summarizing results... \n")
source("scripts/autologistic_04_summary.R")
