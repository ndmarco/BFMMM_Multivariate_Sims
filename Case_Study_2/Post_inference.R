library(ggplot2)
library(ggtern)
library(BayesFMMM)

setwd()

### Load Data
dat <- readRDS("./BRCA_dat_complete.csv")

dir <- paste0(getwd(),"/trace")

### Plot membership
Z <- ZCI(dir, 50, rescale = F)


df = data.frame(x = Z$CI_50[,1],
                y = Z$CI_50[,2],
                z = Z$CI_50[,3],
                Subtype = dat$subtype)

ggtern(data = df, aes(x,y,z, color = Subtype)) + geom_point() + tern_limits(1.02,1.02,1.02) + theme_bw() + xlab("Feature 1") + ylab("Feature 2") + zlab("Feature 3")


### Plot Mean
est_mean <- MVMeanCI(dir, 50, rescale = F)
df = data.frame(Expression = c(est_mean$CI_50[1,], est_mean$CI_50[2,],est_mean$CI_50[3,]),
                y = c(rep("Feature 1", 50), rep("Feature 2", 50), rep("Feature 3", 50)),
                gene = rep(colnames(dat)[3:52], 3))
ggplot(data = df, mapping = aes(x = gene, y = y, fill = Expression)) +
  geom_tile()  + scale_fill_gradient(
    low = "lightblue",
    high = "lightcoral"
  ) + # Add a nicer x-axis title
  theme(axis.title.y = element_blank(), # Remove the y-axis title
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

PAM50 <- read.csv("./pam50_centroids.csv")

df <- data.frame(Expression = c(PAM50$LumA, PAM50$Basal, PAM50$Her2),
                 y = c(rep("LumA", 50), rep("Basal", 50), rep("Her2", 50)),
                 gene = rep(colnames(dat)[3:52], 3))
df$y <- factor(df$y, ordered = T, levels = c("LumA", "Basal", "Her2"))
df$gene <- as.factor(df$gene)
ggplot(data = df, mapping = aes(x = gene, y = y, fill = Expression)) +
  geom_tile()  + scale_fill_gradient(
    low = "lightblue",
    high = "lightcoral"
  ) + # Add a nicer x-axis title
  theme(axis.title.y = element_blank(), # Remove the y-axis title
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

#### Correlation Plots
cov_est1 <- MVCovCI(dir, 50, 1, 1, rescale = F)
corr1 <- cov_est1$CI_50
for(i in 1:nrow(corr1)){
  for(j in 1:ncol(corr1)){
    if( i != j){
    corr1[i,j] <- corr1[i,j] /(sqrt(corr1[i,i] * corr1[j,j]))
    }
    if(i > j){
      corr1[i,j] <- 0
    }
  }
}

colnames(corr1) <- colnames(dat)[3:52]
rownames(corr1) <- colnames(dat)[3:52]

mat <- corr1
df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df <- df[df$from != df$to,]
df <- df[abs(df$value) > 0.8,]
library(circlize)

cols = colorRamp2(c(-1,0,1),c("blue","white","red"),transparency = 0)
chordDiagram(df,col=cols, annotationTrack = c("name","grid"), transparency = 0.2)

cov_est1 <- MVCovCI(dir, 50, 1000, 2, 2, rescale = F)
corr1 <- cov_est1$CI_50
for(i in 1:nrow(corr1)){
  for(j in 1:ncol(corr1)){
    if( i != j){
      corr1[i,j] <- corr1[i,j] /(sqrt(corr1[i,i] * corr1[j,j]))
    }
    if(i > j){
      corr1[i,j] <- 0
    }
  }
}

colnames(corr1) <- colnames(dat)[3:52]
rownames(corr1) <- colnames(dat)[3:52]

mat <- corr1
df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df <- df[df$from != df$to,]
df <- df[abs(df$value) > 0.8,]
library(circlize)

cols = colorRamp2(c(-1,0,1),c("blue","white","red"),transparency = 0)
chordDiagram(df,col=cols, annotationTrack = c("name","grid"), transparency = 0.2)


cov_est1 <- MVCovCI(dir, 50, 3, 3, rescale = F)
corr1 <- cov_est1$CI_50
for(i in 1:nrow(corr1)){
  for(j in 1:ncol(corr1)){
    if( i != j){
      corr1[i,j] <- corr1[i,j] /(sqrt(corr1[i,i] * corr1[j,j]))
    }
    if(i > j){
      corr1[i,j] <- 0
    }
  }
}

colnames(corr1) <- colnames(dat)[3:52]
rownames(corr1) <- colnames(dat)[3:52]

mat <- corr1
df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df <- df[df$from != df$to,]
df <- df[abs(df$value) > 0.8,]
library(circlize)

cols = colorRamp2(c(-1,0,1),c("blue","white","red"),transparency = 0)
chordDiagram(df,col=cols, annotationTrack = c("name","grid"), transparency = 0.2)
