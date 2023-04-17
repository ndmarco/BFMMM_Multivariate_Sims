library(funFEM)
### Real Case study
library(BayesFMMM)
library(mixtools)
library(mclust)

setwd("")

#################################################################
## Change relevant directories and make folders before running ##
#################################################################


### Peak alpha data
library(pracma)
library(gridExtra)
# Subject ID
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
# Channel ID (order of chan_id corresponds to 1:25 labeling of regions)
chan_id <- c('FP1', 'FP2','F9','F7','F3','Fz','F4','F8','F10','T9','T7',
             'C3','CZ','C4','T8','T10','P9','P7','P3','PZ','P4','P8','P10','O1','O2')

# Demographic Data
demDat <- read.csv(file='demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("pa.dat.Rdata")
# ID: subject ID
# group: TD(1) or ASD (2)
# func: frequency domain
# reg: electrode (order corresponds to chan_id above)
# Age: age in months
# y: alpha spectra density
out1 <- unique(pa.dat$func)
out3 <- unique(pa.dat$reg)
matplot(matrix(pa.dat$y, nrow = length(out1)), type = "l") # data
trapz(out1, pa.dat$y[1:33]) # all functional observations integrate to 1 (normalized across electordes, subjects)

### Convert to wide format
Y <- pa.dat
## paper used T8 electrode
Y <- Y[Y$reg == 15,]
Y$ID <- paste(Y$ID, Y$reg, sep = ".")
Y <- reshape(Y[,c(1,3,6)], idvar = "ID", timevar = "func", direction = "wide")
Y <- Y[,-1]
Y <- as.matrix(Y)
matplot(t(Y), type = "l")

##################################
## Fitting a mixture of normals ##
##################################
library(tidyverse)
# function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(Y, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- seq(1,16)

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

model <- kmeans(Y,3)
model$centers
model$cluster
df <- matrix(0, nrow = 99, ncol = 3)
df[1:33,1] <- model$centers[1,]
df[1:33,3] <- 1
df[34:66,1] <- model$centers[2,]
df[34:66,3] <- 2
df[67:99,1] <- model$centers[3,]
df[67:99,3] <- 3
df[,2] <- rep(seq(6,14,0.25), 3)
colnames(df) <- c("Power", "Frequency (Hz)", "Cluster")
df <- as.data.frame(df)
df[,3] <- as.factor(df[,3])
p <-  ggplot(df, aes(x = `Frequency (Hz)`, y = `Power`, color = `Cluster`)) + geom_point() + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab("Relative Power")

df <- matrix(0, nrow = 33*20, ncol = 4)
for(i in 1:20){
  df[((i-1)*33 + 1):(i*33), 1] <- Y[i,]
  df[((i-1)*33 + 1):(i*33), 2] <- seq(6,14,0.25)
  df[((i-1)*33 + 1):(i*33), 3] <- model$cluster[i]
  df[((i-1)*33 + 1):(i*33), 4] <- i
}
colnames(df) <- c("Power", "Frequency (Hz)", "Cluster", "id")
df <- as.data.frame(df)
df[,3] <- as.factor(df[,3])
p1 <-  ggplot(df, aes(x = `Frequency (Hz)`, y = `Power`, color = `Cluster`, group = `id`)) + geom_line() + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab("Relative Power")

grid.arrange(p, p1, nrow = 1)
