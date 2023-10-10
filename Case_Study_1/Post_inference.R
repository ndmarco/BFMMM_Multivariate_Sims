library(BayesFMMM)
library(eegkit)
library(gridExtra)
library(grDevices)
library(ggplot2)
library(MASS)

setwd("")
########################################
## Note: Data not publicly available ###
########################################

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

load("./pa.dat.Rdata")
chan_id <- c('FP1', 'FP2','F9','F7','F3','FZ','F4','F8','F10','T9','T7',
             'C3','CZ','C4','T8','T10','P9','P7','P3','PZ','P4','P8','P10','O1','O2')
chan_id_sub <- c('F5', 'F6', 'T7', 'CZ', 'T8', 'PZ')

setwd("")
dir <- paste0(getwd(), "/")
mean_est <- MVMeanCI(dir, 50)

#### Correlation Plots
library(reshape2)
cov_est1 <- MVCovCI(dir, 50, 1, 1, burnin_prop = 0.5)
mat <- cov_est1$CI_50
rownames(mat) <- colnames(mat) <- seq(6,14, 0.25)
melted_mat <- melt(mat)
ggplot(data = melted_mat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + ggtitle("Covariance for Feature 1") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none") + xlab("Frequency (Hz)") + ylab("Frequency (Hz)")

cov_est2 <- MVCovCI(dir, 50, 2, 2, burnin_prop = 0.5)
mat <- cov_est2$CI_50
rownames(mat) <- colnames(mat) <- seq(6,14, 0.25)
melted_mat <- melt(mat)
ggplot(data = melted_mat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + ggtitle("Covariance for Feature 2") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none") + xlab("Frequency (Hz)") + ylab("Frequency (Hz)")

cov_est12 <- MVCovCI(dir, 50, 1, 2, burnin_prop = 0.5)
mat <- cov_est12$CI_50
rownames(mat) <- colnames(mat) <- seq(6,14, 0.25)
melted_mat <- melt(mat)
ggplot(data = melted_mat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + ggtitle("Cross-Covariance") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none") + xlab("Frequency (Hz)") + ylab("Frequency (Hz)")


#### Mean plots
library(tibble)
dir <- paste0(getwd(), "/")
mean_est <- MVMeanCI(dir, 50, burnin_prop = 0.5)
mean <- c(mean_est$CI_50[1,],mean_est$CI_Upper[1,], mean_est$CI_Lower[1,])
shape <- as.factor(c(rep(1, 33), rep(2,66)))
df_1 <- data_frame("power" = mean, "freq" = rep(seq(6,14,0.25), 3), "shape" = shape)


p1 <- ggplot(data = df_1, aes(x = freq, y = power, shape = shape, color = shape, size = shape)) + geom_point() + scale_shape_manual(values=c(20, 95))  +
  scale_size_manual(values = c(2,4)) + scale_color_manual(values = c("black", "blue")) + ggtitle("Feature 1") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none") + xlab("Frequency (Hz)") + ylab("Relative Power")

mean2 <- c(mean_est$CI_50[2,],mean_est$CI_Upper[2,], mean_est$CI_Lower[2,])
df_2 <- data_frame("power" = mean2, "freq" = rep(seq(6,14,0.25), 3), "shape" = shape)

p2 <- ggplot(data = df_2, aes(x = freq, y = power, shape = shape, color = shape, size = shape)) + geom_point() + scale_shape_manual(values=c(20, 95))  +
  scale_size_manual(values = c(2,4)) + scale_color_manual(values = c("black", "blue")) + ggtitle("Feature 2") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none") + xlab("Frequency (Hz)") + ylab("Relative Power")

grid.arrange(p1, p2, ncol = 2)


### Membership Plot
library(ggplot2)
Z <- ZCI(dir, 50)
demDat <- read.csv(file='./demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]
data_Z <- data.frame("Cluster 1" = Z$CI_50[,1], "Clinical Diagnosis" = demDat$Group)
data_Z$Clinical.Diagnosis[data_Z$Clinical.Diagnosis == 2] <- "ASD"
data_Z$Clinical.Diagnosis[data_Z$Clinical.Diagnosis == 1] <- "TD"
p3 <- ggplot(data= data_Z, aes(x = `Cluster.1` , y = Clinical.Diagnosis)) + geom_violin(trim = F) + geom_point() + xlab("Feature 1 Membership") + ylab("Clinical Diagnosis") +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+ xlim(c(0,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        panel.background = element_blank(),axis.line = element_line(colour = "black"),
                                        plot.title = element_text(hjust = 0.5))

### Plot Trajectory of Means
means <- matrix(0, nrow = 101, ncol= 33)
for(i in 0:100){
  means[i+1,] <- ((i * 0.01) * mean_est$CI_50[1,]) + ((1 - (i * 0.01)) * mean_est$CI_50[2,])
}

df <- data_frame("power" = as.vector(t(means)), "freq" = rep(seq(6,14,0.25), 101), "Z" = rep(0,3333))
for(i in 2:101){
  df$Z[((i-1)*33 + 1):(i*33)] <-((i) * 0.01)
}
library(latex2exp)

p <- ggplot(data = df, aes(x = freq, y = power, color = Z, group = Z)) + geom_point() +
  ggtitle("Mean Structure") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)) + xlab("Frequency (Hz)") + ylab("Relative Power") +
  scale_colour_gradientn(name = TeX("$Z_{i1}$"), colours = c('blue', 'red'))


grid.arrange(p1, p2, p, p3,  layout_matrix = rbind(c(1,2),c(1,2), c(1,2), c(3,4),c(3,4),c(3,4)))

####### Create Plot of Uncertainty

df_1 <- matrix(0, nrow = 10000, ncol = 33)

for(i in 1:10000){
  df_1[i,] <- mvrnorm(n=1, mean_est$mean_trace[1,,14999 +i], cov_est1$cov_trace[,,1499+i])
}

get_quantiles <- function(df){
  df_1_median <- rep(0, 33)
  df_1_upper <- rep(0, 33)
  df_1_lower <- rep(0, 33)

  for(i in 1:33){
    quant <- quantile(df[,i], c(0.25, 0.5 , 0.75))
    df_1_lower[i] <- quant[1]
    df_1_median[i] <- quant[2]
    df_1_upper[i] <- quant[3]
  }
  return(list("lower" = df_1_lower, "upper" = df_1_upper, "median" = df_1_median))
}
quant_1 <- get_quantiles(df_1)

df_2 <- matrix(0, nrow = 10000, ncol = 33)
for(i in 1:10000){
  df_2[i,] <- mvrnorm(n=1, mean_est$mean_trace[2,,14999 +i], cov_est2$cov_trace[,,1499+i])
}
quant_2 <- get_quantiles(df_2)
df_12 <- matrix(0, nrow = 10000, ncol = 33)
for(i in 1:10000){
  df_12[i,] <- mvrnorm(n=1, (mean_est$mean_trace[2,,14999 +i] + mean_est$mean_trace[1,,14999 +i])/2, (0.25 *cov_est2$cov_trace[,,1499+i] + 0.25*cov_est1$cov_trace[,,1499+i] +
                                                                                                        0.25* cov_est12$cov_trace[,,1499+i] + t(0.25* cov_est12$cov_trace[,,1499+i])))
}
quant_12 <- get_quantiles(df_12)


mean <- c(quant_1$median, quant_2$median, quant_12$median)
colors <- as.factor(c(rep(1, 33), rep(2,33), rep(3, 33)))
df <- data_frame("power" = mean, "freq" = rep(seq(6,14,0.25), 3), "colors" = colors)
df_1 <- df[df$colors ==1,]
df_2 <- df[df$colors ==2,]
df_3 <- df[df$colors ==3,]

p1 <- ggplot(data = df, aes(x = freq, y = power, color = colors)) + geom_point() +
  geom_ribbon(data = df_1, aes(ymin = quant_1$lower, ymax = quant_1$upper, x = seq(6,14,0.25)), fill = "red", alpha = 0.2) +
  geom_ribbon(data = df_2, aes(ymin = quant_2$lower, ymax = quant_2$upper, x = seq(6,14,0.25)), fill = "blue", alpha = 0.2) +
  geom_ribbon(data = df_3, aes(ymin = quant_12$lower, ymax = quant_12$upper, x = seq(6,14,0.25)), fill = "purple", alpha = 0.2) +
  scale_color_manual(values = c("red", "blue", "purple")) + ggtitle("Estimated Distribution by Feature Membership") + theme_classic() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5), legend.position = "none") + xlab("Frequency (Hz)") + ylab("Relative Power")


cols <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=101))




### AIC BIC DIC

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

load("./ASD_multivariate/pa.dat.Rdata")
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

dir <- "./trace/"

AIC_2 <- MVAIC(dir, 50, Y)
BIC_2 <- MVBIC(dir, 50, Y)
DIC_2 <- MVDIC(dir, 50, Y)

dir <- "./trace_3/"
AIC_3 <- MVAIC(dir, 50, Y)
BIC_3 <- MVBIC(dir, 50, Y)
DIC_3 <- MVDIC(dir, 50, Y)
