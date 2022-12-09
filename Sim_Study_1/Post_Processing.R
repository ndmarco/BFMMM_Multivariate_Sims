library(BayesFMMM)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(pracma)

### Set working dir
setwd()

err_Z <- matrix(0, 50, 3)
err_mean1 <- matrix(0, 50, 3)
err_mean2 <- matrix(0, 50, 3)
err_cov1 <- matrix(0, 50, 3)
err_cov2 <- matrix(0, 50, 3)
err_cov12 <- matrix(0, 50, 3)

norm_mu1 <- matrix(1, nrow = 50, ncol = 3)
norm_mu2 <- matrix(1, nrow = 50, ncol = 3)
norm_C1 <- matrix(1, nrow = 50, ncol = 3)
norm_C2 <- matrix(1, nrow = 50, ncol = 3)
norm_C12 <- matrix(1, nrow = 50, ncol = 3)

for(j in 1:3){
  if(j == 1){
    dir <- "./50_obs/"
  }
  if(j == 2){
    dir <- "./250_obs/"
  }
  if(j == 3){
    dir <- "./1000_obs/"
  }
  for(i in 1:50){
    x <- readRDS(paste(dir, "sim", i, "/truth.RDS", sep=""))
    x$Phi_true <- x$Phi
    x$nu_true <- x$nu
    Z_true <- x$Z
    nu_1_true <-  x$nu[1,]
    nu_2_true <-  x$nu[2,]
    norm_mu1[i,j] <- norm(as.matrix(x$mu_1))
    norm_mu2[i,j] <- norm(as.matrix(x$mu_2))
    cov1_true <- matrix(0, 10, 10)
    cov2_true <- matrix(0, 10, 10)
    cov12_true <- matrix(0, 10, 10)
    for(m in 1:4){
      cov1_true <- cov1_true + t(t(x$Phi[1, ,i])) %*% x$Phi[1, , i]
      cov2_true <- cov2_true + t(t(x$Phi[2, ,i])) %*% x$Phi[2, , i]
      cov12_true <- cov12_true + t(t(x$Phi[1, ,i])) %*% x$Phi[2, , i]
    }
    
    norm_C1[i,j] <- norm(cov1_true)
    norm_C2[i,j] <- norm(cov2_true)
    norm_C12[i,j] <- norm(cov12_true)

    dir_i <- paste(dir, "sim", i, "/", sep = "")
    est_nu <- MVMeanCI(dir_i, 20, burnin_prop =  0.5)
    norm_1 <- norm(est_nu$CI_50[1,] - as.matrix(x$mu_1))
    norm_2 <- norm(est_nu$CI_50[2,] - as.matrix(x$mu_1))
    Z_est <- ZCI(dir_i, 20, burnin_prop = 0.5)
    if(norm_1 < norm_2){
      err_mean1[i,j] <- norm_1
      err_mean2[i,j] <- norm(est_nu$CI_50[2,] - as.matrix(x$mu_2))
      est_cov <- MVCovCI(dir_i, 20, 100, 1, 1, burnin_prop = 0.5)
      err_cov1[i,j] <- norm(est_cov$CI_50 - cov1_true)
      est_cov <- MVCovCI(dir_i, 20, 100, 2, 2, burnin_prop = 0.5)
      err_cov2[i,j] <- norm(est_cov$CI_50 - cov2_true)
      est_cov <- MVCovCI(dir_i, 20, 100, 1, 2, burnin_prop = 0.5)
      err_cov12[i,j] <- norm(est_cov$CI_50 - cov12_true)
      err_Z[i,j] <- sqrt(norm(Z_est$CI_50 - Z_true, "F") / (2 * nrow(Z_est$CI_50)))
    }else{
      err_mean1[i,j] <- norm_2
      err_mean2[i,j] <- norm(est_nu$CI_50[1,] - as.matrix(x$mu_2))
      est_cov <- MVCovCI(dir_i, 20, 100, 2, 2, burnin_prop = 0.5)
      err_cov1[i,j] <- norm(est_cov$CI_50 - cov1_true)
      est_cov <- MVCovCI(dir_i, 20, 100, 1, 1, burnin_prop = 0.5)
      err_cov2[i,j] <- norm(est_cov$CI_50 - cov2_true)
      est_cov <- MVCovCI(dir_i, 20, 100, 2, 1, burnin_prop = 0.5)
      err_cov12[i,j] <- norm(est_cov$CI_50 - cov12_true)
      Z_est_i <- Z_est
      Z_est$CI_50[,1] <- Z_est_i$CI_50[,2]
      Z_est$CI_50[,2] <- Z_est_i$CI_50[,1]
      err_Z[i,j] <- sqrt(norm(Z_est$CI_50 - Z_true, "F") / (2 * nrow(Z_est$CI_50)))
    }
    print(i)
    print(j)
    
  }
}


C1_RMSE <- matrix(0, 150, 2)
C1_RMSE[1:50,1] <- (int_err_cov1[,1] / norm_C1[,1])
C1_RMSE[1:50,2] <- 50
C1_RMSE[51:100,1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[51:100,2] <- 250
C1_RMSE[101:150,1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[101:150,2] <- 1000
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("RSE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p1 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, lim = c(0,0.5)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
C2_RMSE <- matrix(0, 150, 2)
C2_RMSE[1:50,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:50,2] <- 50
C2_RMSE[51:100,1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[51:100,2] <- 250
C2_RMSE[101:150,1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[101:150,2] <- 1000
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("RSE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p2 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, lim = c(0,0.5)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
C12_RMSE <- matrix(0, 150, 2)
C12_RMSE[1:50,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:50,2] <- 50
C12_RMSE[51:100,1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[51:100,2] <- 250
C12_RMSE[101:150,1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[101:150,2] <- 1000
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("RSE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p3 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, lim = c(0,0.5)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
mu1_RMSE <- matrix(0, 150, 2)
mu1_RMSE[1:50,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:50,2] <- 50
mu1_RMSE[51:100,1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[51:100,2] <- 250
mu1_RMSE[101:150,1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[101:150,2] <- 1000
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("RSE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p4 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, lim= c(0,.05))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
mu2_RMSE <- matrix(0, 150, 2)
mu2_RMSE[1:50,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:50,2] <- 50
mu2_RMSE[51:100,1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[51:100,2] <- 250
mu2_RMSE[101:150,1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[101:150,2] <- 1000
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("RSE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p5 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, lim=c(0,0.05)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))


Z_RMSE <- matrix(0, 150, 2)
Z_RMSE[1:50,1] <- err_Z[,1]
Z_RMSE[1:50,2] <- 50
Z_RMSE[51:100,1] <- err_Z[,2]
Z_RMSE[51:100,2] <- 250
Z_RMSE[101:150,1] <- err_Z[,3]
Z_RMSE[101:150,2] <- 1000
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p6 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") +
  geom_boxplot() +  theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5))
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                           c(4, 4, 5, 5, 6, 6)))