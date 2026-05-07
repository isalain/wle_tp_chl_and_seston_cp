# wle_tp_chl_and_seston_cp
This repository contains the R code on Bayesian hierarchical models (BHMs) of the total phosphorus (TP)–chlorophyll relationship and seston C:P stoichiometry in western Lake Erie, using GLERL–CIGLR water-quality monitoring data archived at NOAA NCEI. The pipeline is split into four sequential snippets that should be run in order. 


Data source
Raw HABs field-sampling data are downloaded directly from the NOAA NCEI archive within 01_data_preprocessing.R:

2012–2018: NCEI accession 0187718
2019: NCEI accession 0209116
2020–2021: NCEI accession 0254720
2022: NCEI accession 0292222

No local data files are required; the script writes the CSVs into the working directory on first run.
How to run

Set the same working directory at the top of every snippet (e.g. work_dir <- "path/to/tp_chl"). Each script begins with setwd(work_dir).
Run the scripts in order:

r   source("01_data_preprocessing.R")
   source("02_model_fitting.R")
   source("03_parameter_extraction.R")
   source("04_figures.R")

Fitted models are cached as model_month.rds, model_site.rds, and seston_model.rds. Subsequent runs reuse the cached models unless the .rds files are deleted.

Software requirements

R ≥ 4.2
Stan toolchain (rstan / cmdstanr backend for brms)
R packages: brms, posterior, bayesplot, dplyr, tidyverse, lubridate, ggplot2, patchwork, ggrepel, ggExtra, viridis, ellipse
