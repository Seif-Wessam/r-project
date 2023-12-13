setwd("C:/Users/fanny/Desktop/stats R/project")

#load useful libraries
library (ade4)
library(vegan)
library(tidyverse)
library(heatmaply)
library(randomForest)
library(gplots)
library(ggplot2)
library(ggpubr)
library("GGally")
library("tidyr")
library(pastecs)
library(psych)
library(NbClust)
library(dendextend)
library(partykit)

data (aravo)

spe<-aravo[[1]]
env<-aravo[[2]]
traits<-aravo[[3]]

##change of variables
#Form
env$Form <- as.numeric(env$Form)

library(dplyr)
env <- env %>%
  mutate(no_ZoogD = ifelse(ZoogD == "no", 1, 0),
         some_ZoogD = ifelse(ZoogD == "some", 1, 0),
         high_ZoogD = ifelse(ZoogD == "high", 1, 0)) %>%
  select(-ZoogD)


#DATA RESEMBLANCE

##Q-mode dissimilarity and distance measures for species
#pq pas le garder, à modifier pour rendre + lisible
#modification of image() to include row and column axis labels and to construct a heatmap of species
par(mfrow=c(1,1))
image(as.matrix(spe)) #original
image.real <- function(mat) { 
  mat <- t(mat)[,nrow(mat):1]
  image(mat, axes = FALSE, col = hcl.colors(15, palette="viridis"))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat), las=2)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat))
  box() 
}
image.real(as.matrix(spe)) #modified
heatmap.2(as.matrix(spe), Colv = FALSE, Rowv =FALSE, dendrogram="none", trace="none")


#percentage difference (Bray-Curtis) dissimilarity matrix on raw species
spe.db <- vegdist(spe)
image.real(as.matrix(spe.db))


#on garde 2 dissimilarity matrices -> subplot joli, pour montrer clusters et qu'elles sont similaires
par(mfrow=c(1,1))

#percentage difference (Bray-Curtis) dissimilarity matrix on log-transformed abundances
spe.dbln <- vegdist(log1p(spe))
image.real(as.matrix(spe.dbln))

#Hellinger distance matrix - one type of euclidian distance
spe.hel <- decostand(spe, "hel")
spe.dh <- dist(spe.hel)
image.real(as.matrix(spe.dh))

# Unsupervised classification method -----------------------------
# quel est le best linkage? trade off cophenetic correlation, compact + same size groups

# Compute matrix of chord distances among sites
spe.norm <- decostand(spe, "normalize")
spe.ch <- vegdist(spe.norm, "euc")

# Compute and plot complete-linkage agglomerative clustering
spe.ch.complete <- hclust(spe.ch, method = "ward.D2")
par(mfrow = c(2, 2))
par(mar = c(2, 5, 2, 2))
plot(spe.ch.complete, main = "Ward linkage")

# Complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)

# Shepard diagram
plot(spe.ch, spe.ch.comp.coph,
     xlab = "Chord distance",
     ylab = "Cophenetic distance",
     asp = 1, xlim = c(0, 2),
     ylim = c(0, sqrt(2)),
     main = c("Shepard diagram", paste("Cophenetic correlation =", round(cor(spe.ch, spe.ch.comp.coph), 3))))
abline(0, 1)
lines(lowess(spe.ch, spe.ch.comp.coph), col = "red", lwd=3)


# find optimal number of clusters 
library(NbClust)

Nb.complete<-NbClust(spe, diss=spe.ch, distance = NULL, min.nc=3, max.nc=10, 
                     method = "complete", index = "ch")
Nb.complete #we find 4 clusters for an optimal choice
plot(3:10,Nb.complete$All.index, xlab="number of clusters", ylab="Calinski and Harabasz index", main="Optimal number of clusters")
abline(v=4, col="red", lty=2)

#plot the dendrogram with these 4 optimal cluster in evidence
library(dendextend)
complete.dend <- as.dendrogram(spe.ch.complete)
colors_to_use <- Nb.complete$Best.partition
colors_to_use<-colors_to_use[order.dendrogram(complete.dend)]
labels_colors(complete.dend)<-1
complete.dend <- complete.dend %>% color_branches(k = 4)
plot(complete.dend, main="Complete linkage with 4 clusters")

#SUPERVISED CLASSIF
spe.cw.g <- cutree(spe.ch.complete, 4)
spe.cw.g

D<-cbind(spe.cw.g,env,spe)
class.groups=ctree(as.factor(spe.cw.g)~ PhysD+Slope+no_ZoogD+some_ZoogD+high_ZoogD+Aspect+Form+Snow, data = D)
plot(class.groups)

class.snow<-ctree(Snow~PhysD+Slope+high_ZoogD+Aspect+Form, data = D)
plot(class.snow)

class.zoo<-ctree(high_ZoogD~Snow+PhysD+Slope+Aspect+Form, data = D)
plot(class.zoo)

rf <- randomForest(as.factor(spe.cw.g)~., env, ntree=500, mtry=, importance=TRUE,
                   na.action=na.omit,do.trace=100,proximity=T)
plot(rf)

#partial plot
par(mfrow=c(2,3))
partialPlot(rf, env, Aspect)
partialPlot(rf, env, Slope)
partialPlot(rf, env, Form)
partialPlot(rf, env, PhysD)
partialPlot(rf, env, high_ZoogD)
partialPlot(rf, env, Snow)

# classification of traits
traits.norm <- decostand(traits, "normalize")
traits.ch <- vegdist(traits.norm, "euc")
tra.ch.ward <- hclust(traits.ch, method = "ward.D2")
plot(tra.ch.ward,  main = "Chord - Ward")
tra.cw.g <- cutree(tra.ch.ward, 4)

# Compute CA
spe.ca <- cca(spe)
summary(spe.ca)		# default scaling 2
summary(spe.ca, scaling = 1)

# CA biplots
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of species
plot(spe.ca, 
     scaling = 1, 
     main = "CA abundances - biplot scaling 1"
)
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main = "CA abundances - biplot scaling 2")

#to remove the horseshoe we use DCA
spe.DCA<-decorana(spe, iweigh=0, iresc=4, ira=0, mk=26, short=0, before=NULL, after=NULL)
spe.DCA
plot(spe.DCA, disp="sites")

#RDA over env variables
env.pca <-rda(env, scale = TRUE)
#env.pca <-rda(env, scale = TRUE, row.wt = spe.ca$lw) #ADDED THE WEIGHT HERE
env.pca
summary(env.pca)
par(mfrow = c(1, 2))
biplot(env.pca, main = "PCA - scaling 2")
biplot(env.pca, scaling = 1, main = "PCA - scaling 1")

#RDA constrained by Hel. matrix
spe.rda <- rda(spe.hel ~ ., env) #"." means constrained for all the env par in env3
summary(spe.rda)

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

## Triplots of the rda results (lc site scores)
## Site scores as linear combinations of the environmental variables
par(mfrow=c(1,3))
plot(spe.rda, scaling = 1,   display = c("lc"), main = "RDA - sites")
plot(spe.rda, scaling = 1,   display = c("sp"), main = "RDA - species")
plot(spe.rda, scaling = 1,   display = c("cn"), main = "RDA - constraints")

spe.good <- goodness(spe.rda)
sel.sp <- which(spe.good[, 2] >= 0.6)
anova(spe.rda, permutations = how(nperm = 999))
vif.cca(spe.rda)

mod0 <- rda(spe.hel ~ 1, data = env)
spe.rda.all <- rda(spe.hel ~ ., data = env)
step.forward <- ordistep(mod0, scope = formula(spe.rda.all), direction = "forward", permutations = how(nperm = 999))
RsquareAdj(step.forward)

traits.pca <- rda(traits, scale = TRUE)
#traits.pca <- rda(traits, scale = TRUE, row.wt = spe.ca$cw) #ADDED THE WEIGHT HERE
traits.pca
summary(traits.pca)
par(mfrow = c(1, 2))

biplot(traits.pca, main = "PCA - scaling 2")
biplot(traits.pca, scaling = 1, main = "PCA - scaling 1")
arrows(0, 0,  traits.pca[, 1], traits.pca[, 2])
text(traits.pca, disp = 'species', scaling = 1, col = "red")

#plot pca env and traits together
par(mfrow = c(1, 2))
biplot(env.pca, scaling = 1, main = "Env PCA")
biplot(traits.pca, scaling = 1, main = "Traits PCA")
arrows(0, 0,  traits.pca[, 1], traits.pca[, 2])
text(traits.pca, disp = 'species', scaling = 1, col = "red")

#try nmds or to filter before cca!!
spe.nmds <- metaMDS(spe, distance = "bray")
spe.nmds$stress
plot(spe.nmds, type = "t",main = paste("NMDS Bray Curtis; Stress =",round(spe.nmds$stress, 3)))

par(mfrow = c(1, 2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
gof
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof * 300) #to see the goodness of fit if the sites




# Unsupervised classification method - test with traits----------------------------


# Compute matrix of chord distances among sites
traits.norm <- decostand(traits, "normalize")
traits.ch <- vegdist(traits.norm, "euc")

# Compute and plot complete-linkage agglomerative clustering
traits.ch.ward <- hclust(traits.ch, method = "ward.D2")
par(mfrow = c(1, 1))
#par(mar = c(2, 5, 2, 2))
plot(traits.ch.ward, main = "Ward linkage")

# ward linkage clustering
traits.ch.comp.coph <- cophenetic(traits.ch.ward)
cor(traits.ch, traits.ch.comp.coph)

# Shepard diagram
plot(traits.ch, traits.ch.comp.coph,
     xlab = "Chord distance",
     ylab = "Cophenetic distance",
     asp = 1, xlim = c(0, 2),
     ylim = c(0, sqrt(2)),
     main = c("Shepard diagram", paste("Cophenetic correlation =", round(cor(spe.ch, spe.ch.comp.coph), 3))))
abline(0, 1)
lines(lowess(spe.ch, spe.ch.comp.coph), col = "red", lwd=3)


# find optimal number of clusters 
# best = 2, then 3, 4 and 9
library(NbClust)

Nb.complete<-NbClust(traits, diss=traits.ch, distance = NULL, min.nc=3, max.nc=20, 
                     method = "complete", index = "ch")
Nb.complete 
plot(3:20,Nb.complete$All.index, xlab="number of clusters", ylab="Calinski and Harabasz index", main="Optimal number of clusters")
#abline(v=2, col="red", lty=2)

#plot the dendrogram with these 4 optimal cluster in evidence
library(dendextend)
complete.dend <- as.dendrogram(traits.ch.ward)
colors_to_use <- Nb.complete$Best.partition
colors_to_use<-colors_to_use[order.dendrogram(complete.dend)]
labels_colors(complete.dend)<-1
complete.dend <- complete.dend %>% color_branches(k = 9)
plot(complete.dend, main="Complete linkage with 4 clusters")

#SUPERVISED CLASSIF
traits.cw.g <- cutree(traits.ch.ward, 3)
traits.cw.g
spe.transpose <- t(spe)
D<-cbind(traits.cw.g, spe.transpose, traits)
#class.groups=ctree(as.factor(traits.cw.g)~ Height+Spread+Angle+Area+Thick+SLA+N_mass+Seed, data = D)
class.groups=ctree(as.factor(traits.cw.g)~ AR07+AR70, data = D)
plot(class.groups)

env.clusters.snow <- env
env.clusters.snow$early_melting <- ifelse(env.clusters.snow$Snow <= 150, 1, 0)
env.clusters.snow$late_melting <- ifelse(env.clusters.snow$Snow > 180, 1, 0)
env.clusters.snow$middle_melting <- ifelse(env.clusters.snow$Snow <= 180 & env.clusters.snow$Snow > 150 , 1, 0)
env.clusters.snow <- subset(env.clusters.snow, select = -Snow)

#RDA over env variables
env.pca <-rda(env.clusters.snow, scale = TRUE)
#env.pca <-rda(env, scale = TRUE, row.wt = spe.ca$lw) #ADDED THE WEIGHT HERE
env.pca
summary(env.pca)
par(mfrow = c(1, 2))
biplot(env.pca, main = "PCA - scaling 2")
biplot(env.pca, scaling = 1, main = "PCA - scaling 1")

# Nonmetric multidimensional scaling (NMDS)
spe.nmds <- metaMDS(spe, distance = "bray")
#if best solution not repeated -> good idea to run it again
spe.nmds
spe.nmds$stress
plot(spe.nmds, type = "t",main = paste("NMDS Bray Curtis; Stress =",round(spe.nmds$stress, 3)))
