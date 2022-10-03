# Warped Gradient-Enhanced Gaussian Process Surrogate Models for Inference with Intractable Likelihoods

This folder provides reproducible code for the manuscript titled *Warped Gradient-Enhanced Gaussian Process Surrogate Models for Inference with Intractable Likelihoods* by Quan Vu, Matthew T. Moores, and Andrew Zammit-Mangion.

## Folders

The folder `data/` stores the NVDI data used in the manuscript.

The folder `deepsurrogate/` stores the code to fit and predict with the surrogate models.

The folder `figures/` stores the figures produced in the manuscript.

The folder `results/` stores the results produced in the manuscript.

The folder `scripts/` stores the code to reproduce the results in the manuscript.

## Instructions

To reproduce the results and the figures in the manuscript:

Ensure the required R packages are installed. The required packages for these scripts are `bayesImageS`, `coda`, `devtools`, `dplyr`, `doParallel`, `ggplot2`, `GiRaF`, `gridExtra`, `PAWL`, and `tensorflow`. The `R` package `PAWL` can be downloaded from https://cran.r-project.org/src/contrib/Archive/PAWL/. All the other packages can be downloaded using `install.packages()`.

A separate .md file contains the instruction to install TensorFlow v1 (that is required to run the scripts) for Linux Debian 64-bit. If TensorFlow is installed in a conda environment, load this environment before running the scripts.

Set `Vu_Moores_ZammitMangion_2021/` as the working directory. Run through the script files in the folder `scripts/` to reproduce either numerical results or figures in the manuscript.

To reproduce the results in Section 4.1 (the example with the Potts model), run the script:

`run_potts.R`

To reproduce the results in Section 4.2 (the example with the hidden Potts model), run the script:

`run_hidden_potts.R`

To reproduce the results in Section 4.3 (the example with the autologistic model), run the script:

`run_autologistic.R`
