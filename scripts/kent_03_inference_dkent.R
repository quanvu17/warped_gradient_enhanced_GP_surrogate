library("dplyr")
library("ggplot2")
library("Directional")
library("gridExtra")

## Clear memory
rm(list=ls())
load("results/kent_simulated_data.rda")

kappa_true <- 2
beta_true <- 0.5

gamma_1 <- c(1, 0, 0)
gamma_2 <- c(0, 1, 0)
gamma_3 <- c(0, 0, 1)
G <- cbind(gamma_1, gamma_2, gamma_3)

## Kent distribution: Set up grids
dkappa <- 0.05
dbeta <- 0.01

kappa_fine <- seq(0.2, 10, by = dkappa)
beta_fine <- seq(0.1, 5, by = dbeta)

kappa_beta_grid <- 
  expand.grid(kappa = kappa_fine, beta = beta_fine) %>%
  filter(beta < kappa/2) %>%
  mutate(logden = NA)

## Evaluate log-density
for(i in 1:nrow(kappa_beta_grid)) {
  kappa_beta_grid$logden[i] <-  dkent(y = y, 
                                   G = G,
                                   param = c(kappa_beta_grid$kappa[i],
                                             kappa_beta_grid$beta[i]),
                                   logden = TRUE) %>% sum()
}

## Exp-normalise density
exp_normalise <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}
kappa_beta_grid$dens <- exp_normalise(kappa_beta_grid$logden) / (dbeta * dkappa)

## Find marginals
kappa_posterior_unnorm <- group_by(kappa_beta_grid, kappa) %>%
                          summarise(dens_unnorm = sum(dens) * dbeta)  %>% 
                          mutate(dens = dens_unnorm / sum(dens_unnorm))

beta_posterior_unnorm <- group_by(kappa_beta_grid, beta) %>%
                         summarise(dens_unnorm = sum(dens)*dkappa) %>% 
                         mutate(dens = dens_unnorm / sum(dens_unnorm))
