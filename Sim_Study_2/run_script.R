library(MASS)
library(DirichletReg)
library(BayesFPMM)

for(q in 1:10){
  n_obs = 200
  nu <- matrix(0,nrow =3, ncol = 20)
  p <- diag(1, 20)

  nu[1,] <- mvrnorm(n=1, mu = rep(0, 20), Sigma = 10*p)
  nu[2,] <- mvrnorm(n=1, mu = rep(0, 20), Sigma = 10*p)
  nu[3,] <- mvrnorm(n=1, mu = rep(0, 20), Sigma = 10*p)
  Phi_1 <- rnorm(60, 0, 1)
  Phi_2 <- rnorm(60, 0, 0.5)
  Phi_3 <- rnorm(60, 0, 0.2)
  Phi <- array(0, dim = c(3, 20, 3))
  Phi[,,1] <- matrix(Phi_1, nrow = 3)
  Phi[,,2] <- matrix(Phi_2, nrow = 3)
  Phi[,,3] <- matrix(Phi_3, nrow = 3)


  chi <- matrix(rnorm(n_obs *3, 0, 1), ncol = 3, nrow=n_obs)

  Z <- matrix(0, nrow = n_obs, ncol = 3)
  alpha <- c(30, 1, 1)
  for(i in 1:(n_obs * 0.2)){
    Z[i,] <- rdirichlet(1, alpha)
  }
  alpha <- c(1, 30, 1)
  for(i in (n_obs * 0.2 + 1):(n_obs * 0.4)){
    Z[i,] <- rdirichlet(1, alpha)
  }
  alpha <- c(1, 1, 30)
  for(i in (n_obs * 0.4 + 1):(n_obs * 0.6)){
    Z[i,] <- rdirichlet(1, alpha)
  }
  alpha <- c(1, 1, 1)
  for(i in (n_obs * 0.6 + 1):n_obs){
    Z[i,] <- rdirichlet(1, alpha)
  }

  y <- matrix(0, nrow = n_obs, ncol = 20)
  for(i in 1:n_obs){
    mean = rep(0,10)
    for(j in 1:3){
      mean = mean + Z[i,j] * nu[j,]
      for(m in 1:3){
        mean = mean + Z[i,j] * chi[i,m] * Phi[j, ,m]
      }
    }
    y[i,] = mvrnorm(n = 1, mean, diag(0.01, 20))
  }

  x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi)
  saveRDS(x, paste("/Users/nicholasmarco/Box Sync/BayesFPMM_Supporting_Files/BMPMM_simulation/Optimal_K/data/data", q, ".RDS", sep = ""))
  Y <- y


  ####### 2 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 2
  n_eigen <- 3
  dir <- "/Users/nicholasmarco/Box Sync/BayesFPMM_Supporting_Files/BMPMM_simulation/Optimal_K/2_clusters/"

  ## Get Estimates of Z and nu
  est1 <- BMVPMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BMVPMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen, est1$Z, est1$nu)
  dir_i <- paste(dir, "trace", q, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BMVPMM_warm_start(tot_mcmc_iters, k, Y, n_eigen,
                                 est1$Z, est1$pi, est1$alpha_3,
                                 est2$delta, est2$gamma, est2$Phi, est2$A,
                                 est1$nu, est1$tau, est2$sigma, est2$chi, dir = dir_i,
                                 thinning_num = 10, r_stored_iters = 10000)

  ####### 3 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 3
  n_eigen <- 3
  dir <- "/Users/nicholasmarco/Box Sync/BayesFPMM_Supporting_Files/BMPMM_simulation/Optimal_K/3_clusters/"

  ## Get Estimates of Z and nu
  est1 <- BMVPMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BMVPMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen, est1$Z, est1$nu)
  dir_i <- paste(dir, "trace", q, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BMVPMM_warm_start(tot_mcmc_iters, k, Y, n_eigen,
                                 est1$Z, est1$pi, est1$alpha_3,
                                 est2$delta, est2$gamma, est2$Phi, est2$A,
                                 est1$nu, est1$tau, est2$sigma, est2$chi, dir = dir_i,
                                 thinning_num = 10, r_stored_iters = 10000)

  ####### 4 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 4
  n_eigen <- 3
  dir <- "/Users/nicholasmarco/Box Sync/BayesFPMM_Supporting_Files/BMPMM_simulation/Optimal_K/4_clusters/"

  ## Get Estimates of Z and nu
  est1 <- BMVPMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BMVPMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen, est1$Z, est1$nu)
  dir_i <- paste(dir, "trace", q, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BMVPMM_warm_start(tot_mcmc_iters, k, Y, n_eigen,
                                 est1$Z, est1$pi, est1$alpha_3,
                                 est2$delta, est2$gamma, est2$Phi, est2$A,
                                 est1$nu, est1$tau, est2$sigma, est2$chi, dir = dir_i,
                                 thinning_num = 10, r_stored_iters = 10000)

  ####### 5 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 5
  n_eigen <- 3
  dir <- "/Users/nicholasmarco/Box Sync/BayesFPMM_Supporting_Files/BMPMM_simulation/Optimal_K/5_clusters/"

  ## Get Estimates of Z and nu
  est1 <- BMVPMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, n_eigen)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BMVPMM_Theta_est(tot_mcmc_iters, n_try, k, Y, n_eigen, est1$Z, est1$nu)
  dir_i <- paste(dir, "trace", q, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BMVPMM_warm_start(tot_mcmc_iters, k, Y, n_eigen,
                                 est1$Z, est1$pi, est1$alpha_3,
                                 est2$delta, est2$gamma, est2$Phi, est2$A,
                                 est1$nu, est1$tau, est2$sigma, est2$chi, dir = dir_i,
                                 thinning_num = 10, r_stored_iters = 10000)
}

