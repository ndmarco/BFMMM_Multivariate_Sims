library(BayesFMMM)

## Before running, set working directory to file location
setwd()

### Load Data
dat <- readRDS("~/BRCA_dat_complete.csv")
Y <- dat[, 3:52]
Y <- as.matrix(Y)
## Set Hyperparameters
tot_mcmc_iters <- 10000
n_try <- 50
k <- 3
n_eigen <- 4

## Run function
est1 <- BMVMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)

n_try <- 6
## Run function
est2 <- BMVMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen, est1$Z, est1$nu)

tot_mcmc_iters <- 500000
dir.create(paste0(getwd(),"/trace"))
MCMC.chain <-BMVMMM_warm_start(tot_mcmc_iters, 3, Y, n_eigen,
                               est1$Z, est1$pi, est1$alpha_3,
                               est2$delta, est2$gamma, est2$Phi, est2$A,
                               est1$nu, est1$tau, est2$sigma, est2$chi, dir = paste0(getwd(),"/trace/"),
                               thinning_num = 10, r_stored_iters = 10000)
