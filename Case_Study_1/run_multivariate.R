library(BayesFPMM)
library(eegkit)

setwd("/Users/user/Box Sync/BayesFMMM_Supporting_Files/ASD_multivariate")

#################################################################
## Change relevant directories and make folders before running ##
#################################################################

subj_id <- sort(c(10,	11,	13,	14,	15,	23,	26,	30,	31,	35,	48,	49,	50,
                  53,	54,	55,	161,165,	184,	188,	189,	195,	201,
                  # 202,	excluded due to low counts
                  207,	210,	213,	214,	242,	255,	261,	282,	283,
                  284,	286,	287,	289,	290,	343,	351,	2,	3,	5,	6,
                  7,	8,	9,	12,	18,	19,	22,	24,	25,	27,	33,	34,	37,	38,
                  40,	41,	42,	43,	44,	47,	51,	401,	405,	406,	408,	411,
                  415,	416,	417,	418,	423,	426,	427,	430,
                  #431,	excluded due to low counts
                  433,	436,	438,	439,	440,	442,	444,	445,	446,	447,
                  448,	450,	451,	452,	453,	3019,	3024,	3026,	3029,	3032))

load("pa.dat.Rdata")
chan_id <- c('FP1', 'FP2','F9','F7','F3','FZ','F4','F8','F10','T9','T7',
             'C3','CZ','C4','T8','T10','P9','P7','P3','PZ','P4','P8','P10','O1','O2')

data <- matrix(0,97,33)
freq <- seq(6,14,0.25)
for(i in 1:length(subj_id)){
  dat_i <- pa.dat[pa.dat$ID == subj_id[[i]], ]
  for(j in 1:length(freq)){
    data_ij <- dat_i[dat_i$func == freq[j],]
    data[i,j] = mean(data_ij$y)
  }
}

Y <- data

## Set Hyperparameters
tot_mcmc_iters <- 8000
n_try <- 50
k <- 2
n_eigen <- 5

## Run function
est1 <- BMVPMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)

## Run function
est2 <- BMVPMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen, est1$Z, est1$nu)

tot_mcmc_iters <- 500000

dir = "/Users/user/Box Sync/BayesFPMM_Supporting_Files/ASD_multivariate/Multivariate/trace/"

MCMC.chain <-BMVPMM_warm_start(tot_mcmc_iters, k, Y, n_eigen,
                               est1$Z, est1$pi, est1$alpha_3,
                               est2$delta, est2$gamma, est2$Phi, est2$A,
                               est1$nu, est1$tau, est2$sigma, est2$chi, thinning_num = 10, r_stored_iters = 10000,
                               dir = dir)

