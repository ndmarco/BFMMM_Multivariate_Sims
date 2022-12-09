library(BayesFMMM)
library(MASS)
library(DirichletReg)
library(future.apply)

run_sim_mean_adj <- function(iter){
  set.seed(iter)
  n_obs_vec <- c(50, 250, 1000)
  for(n in 1:3){
    dir.create(paste0(n_obs_vec[n],"_obs"))
    dir.create(paste0(n_obs_vec[n],"_obs/sim",iter))
    n_obs <- n_obs_vec[n]
    
    nu <- rnorm(20, 0, 3)
    nu <- matrix(nu, nrow = 2, ncol = 10)
    decomp <- svd(nu, nv = 10)
    Phi_i <- matrix(0, nrow = 2, ncol = 10)
    Phi_i[1,] <- t(rnorm(8, 0, 1)) %*% t(decomp$v[,3:10])
    Phi_i[2,] <- t(rnorm(8, 0, 1)) %*% t(decomp$v[,3:10])
    Phi_1 <- Phi_i
    Phi_i[1,] <- t(rnorm(8, 0, 0.7)) %*% t(decomp$v[,3:10])
    Phi_i[2,] <- t(rnorm(8, 0, 0.7)) %*% t(decomp$v[,3:10])
    Phi_2 <- Phi_i
    Phi_i[1,] <- t(rnorm(8, 0, 0.5)) %*% t(decomp$v[,3:10])
    Phi_i[2,] <- t(rnorm(8, 0, 0.5)) %*% t(decomp$v[,3:10])
    Phi_3 <- Phi_i
    Phi_i[1,] <- t(rnorm(8, 0, 0.3)) %*% t(decomp$v[,3:10])
    Phi_i[2,] <- t(rnorm(8, 0, 0.3)) %*% t(decomp$v[,3:10])
    Phi_4 <- Phi_i
    Phi <- array(0, dim = c(2, 10, 4))
    Phi[,,1] <- matrix(Phi_1, nrow = 2)
    Phi[,,2] <- matrix(Phi_2, nrow = 2)
    Phi[,,3] <- matrix(Phi_3, nrow = 2)
    Phi[,,4] <- matrix(Phi_4, nrow = 2)
    
    chi <- matrix(rnorm(n_obs *4, 0, 1), ncol = 4, nrow=n_obs)
    
    Z <- matrix(0, nrow = n_obs, ncol = 2)
    alpha <- c(10, 1)
    for(i in 1:(n_obs * 0.3)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 10)
    for(i in (n_obs * 0.3 + 1):(n_obs * 0.6)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 1)
    for(i in (n_obs * 0.6 + 1):n_obs){
      Z[i,] <- rdirichlet(1, alpha)
    }
    
    
    y <- matrix(0, nrow = n_obs, ncol = 10)
    for(i in 1:n_obs){
      mean = rep(0,10)
      for(j in 1:2){
        mean = mean + Z[i,j] * nu[j,]
        for(m in 1:4){
          mean = mean + Z[i,j] * chi[i,m] * Phi[j, ,m]
        }
      }
      y[i,] = mvrnorm(n = 1, mean, diag(0.01, 10))
    }
    
    x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi)
    
    saveRDS(x, paste("./", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))
    
    
    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 50
    k <- 2
    n_eigen <- 4
    Y <-y
    
    ## Get Estimates of Z and nu
    est1 <- BMVMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)
    
    tot_mcmc_iters <- 4000
    n_try <- 5
    ## Run function
    est2 <- BMVMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen,
                             est1$Z, est1$nu)
    dir_i <- paste("./", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 200000
    MCMC.chain <- BMVMMM_warm_start(tot_mcmc_iters, k, Y, n_eigen, est1$Z, est1$pi, est1$alpha_3,
                                                est2$delta, est2$gamma, est2$Phi, est2$A,
                                                est1$nu, est1$tau, est2$sigma, est2$chi,
                                                dir = dir_i, thinning_num = 100, r_stored_iters = 10000)
  }
}



##### Run Simulation

### Set working dir
setwd("/Users/nicholasmarco/Box Sync/BayesFPMM_Supporting_Files/BMPMM_simulation/")

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)

already_ran <- dir(paste0(getwd(), "/1000_obs"))
to_run <- which(!paste0("sim", 1:25) %in% already_ran)
seeds <- to_run
future_lapply(seeds, function(this_seed) run_sim_mean_adj(this_seed))
