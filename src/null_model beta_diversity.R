# Ecological diversity analyses codes originally written by Meng Chen during the doctoral dissertation
# 05Oct2016
# revised on 12Jan2016
# reorganized on 23Oct2019

###########################################################
# Null model testing and beta diversity calculation       #
# Author: Meng Chen    Date: 05Oct2015   Updated12Jan2016 #
###########################################################

### libraries used in the data analyses
library(permute)
library(lattice)
library(vegan)
library(MASS)
library(EcoSimR)

# set up directory
setwd("your directory")

#################################################
# Parameters and Functions used in the analyses #
#################################################

# occurrence data function
eco.occurrence.func <- function(eco.data) {
  eco.occurrence <- eco.data[,c(1,3,4,5)]
  eco.occurrence$combination <- paste(eco.occurrence[,2], eco.occurrence[,3], eco.occurrence[,4])
  eco.occurrence <- eco.occurrence[,c(1,5)]
  eco.occurrence <- t(table(eco.occurrence))
  eco.occurrence[eco.occurrence > 0 ] <- 1
  return(eco.occurrence)
}

# null model testing function
cooc.null.func <- function (input.data) {
  output.data <- cooc_null_model(input.data, algo = "sim9", metric = "c_score", nReps = 1000, burn_in = 500)
  return(output.data)
}

######################################################
# Null model testing based on different environments #
######################################################

# Read the data of small-bodied mammal communities
mcm.habitat <- read.csv(file="small mammal community data", sep=",", header=T)

# All communities
# ecological parameter combination occurences data
eco.occurrence.all <- eco.occurrence.func(mcm.habitat)
cooc.null.all <- cooc.null.func(eco.occurrence.all)
summary(cooc.null.all)
plot(cooc.null.all, type="hist")

# 2 habitats
# data
mcm.open <- mcm.habitat[mcm.habitat$Habitat == "Open",]
mcm.close <- mcm.habitat[mcm.habitat$Habitat == "Close",]
# converted to occurrence data
eco.occurrence.open <- eco.occurrence.func(mcm.open)
eco.occurrence.close <- eco.occurrence.func(mcm.close)
# calculating C-score
cooc.null.open <- cooc.null.func(eco.occurrence.open)
cooc.null.close <- cooc.null.func(eco.occurrence.close)
# summary
summary(cooc.null.open)
summary(cooc.null.close)
# plots
plot(cooc.null.open, type="hist")
plot(cooc.null.close, type="hist")

# 4 climates
# data
mcm.tropical.climate <- mcm.habitat[mcm.habitat$Climate == "Tropical",]
mcm.arid.climate <- mcm.habitat[mcm.habitat$Climate == "Arid",]
mcm.temperate.climate <- mcm.habitat[mcm.habitat$Climate== "Temperate",]
mcm.cold.climate <- mcm.habitat[mcm.habitat$Climate == "Cold",]
# converted to occurrence data
eco.occurrence.tropical.climate <- eco.occurrence.func(mcm.tropical.climate)
eco.occurrence.arid.climate <- eco.occurrence.func(mcm.arid.climate)
eco.occurrence.temperate.climate <- eco.occurrence.func(mcm.temperate.climate)
eco.occurrence.cold.climate <- eco.occurrence.func(mcm.cold.climate)
# calculating C-score
cooc.null.tropical.climate <- cooc.null.func(eco.occurrence.tropical.climate)
cooc.null.arid.climate <- cooc.null.func(eco.occurrence.arid.climate)
cooc.null.temperate.climate <- cooc.null.func(eco.occurrence.temperate.climate)
cooc.null.cold.climate <- cooc.null.func(eco.occurrence.cold.climate)
# summary
summary(cooc.null.tropical.climate)
summary(cooc.null.arid.climate)
summary(cooc.null.temperate.climate)
summary(cooc.null.cold.climate)
# plots
plot(cooc.null.tropical.climate, type="hist")
plot(cooc.null.arid.climate, type="hist")
plot(cooc.null.temperate.climate, type="hist")
plot(cooc.null.cold.climate, type="hist")

# 8 vegetations
# data
mcm.boreal.forest <- mcm.habitat[mcm.habitat$Vegetation == "Boreal forest",]
mcm.temperate.desert <- mcm.habitat[mcm.habitat$Vegetation == "Desert",]
mcm.temperate.forest <- mcm.habitat[mcm.habitat$Vegetation == "Temperate forest",]
mcm.temperate.grassland <- mcm.habitat[mcm.habitat$Vegetation== "Grassland",]
mcm.shrubland <- mcm.habitat[mcm.habitat$Vegetation == "Shrubland",]
mcm.savanna <- mcm.habitat[mcm.habitat$Vegetation == "Savanna",]
mcm.tropical.seasonal.forest <- mcm.habitat[mcm.habitat$Vegetation == "Tropical seasonal forest",]
mcm.tropical.rainforest <- mcm.habitat[mcm.habitat$Vegetation == "Tropical rainforest",]
# converted to occurrence data
eco.occurrence.tropical.rainforest <- eco.occurrence.func(mcm.tropical.rainforest)
eco.occurrence.tropical.seasonal.forest <- eco.occurrence.func(mcm.tropical.seasonal.forest)
eco.occurrence.savanna <- eco.occurrence.func(mcm.savanna)
eco.occurrence.shrubland <- eco.occurrence.func(mcm.shrubland)
eco.occurrence.desert <- eco.occurrence.func(mcm.temperate.desert)
eco.occurrence.temperate.forest<- eco.occurrence.func(mcm.temperate.forest)
eco.occurrence.boreal.forest <- eco.occurrence.func(mcm.boreal.forest)
eco.occurrence.grassland <- eco.occurrence.func(mcm.temperate.grassland)
# calculating C-score
cooc.null.tropical.rainforest <- cooc.null.func(eco.occurrence.tropical.rainforest)
cooc.null.tropical.seasonal.forest <- cooc.null.func(eco.occurrence.tropical.seasonal.forest)
cooc.null.savanna <- cooc.null.func(eco.occurrence.savanna)
cooc.null.shrubland <- cooc.null.func(eco.occurrence.shrubland)
cooc.null.grassland <- cooc.null.func(eco.occurrence.grassland)
cooc.null.desert <- cooc.null.func(eco.occurrence.desert)
cooc.null.temperate.forest <- cooc.null.func(eco.occurrence.temperate.forest)
cooc.null.boreal.forest <- cooc.null.func(eco.occurrence.boreal.forest)
# summary
summary(cooc.null.tropical.rainforest)
summary(cooc.null.tropical.seasonal.forest)
summary(cooc.null.savanna)
summary(cooc.null.shrubland)
summary(cooc.null.grassland)
summary(cooc.null.desert)
summary(cooc.null.temperate.forest)
summary(cooc.null.boreal.forest)
# plots
plot(cooc.null.tropical.rainforest, type="hist")
plot(cooc.null.tropical.seasonal.forest, type="hist")
plot(cooc.null.savanna, type="hist")
plot(cooc.null.shrubland, type="hist")
plot(cooc.null.grassland, type="hist")
plot(cooc.null.desert, type="hist")
plot(cooc.null.temperate.forest, type="hist")
plot(cooc.null.boreal.forest, type="hist")

########################################################
# Ecological beta diversity (multivariate disperson) #
########################################################

# data
community.matrix <- eco.occurrence
community.matrix <- t(community.matrix)
enviro.data <- unique(mcm.habitat[,c(1:2, 6:7)])

# Based on Oksanen 2011, 2016 approach using jaccard similarity index
# 1 "w" = (b+c)/(2*a+b+c)
# 2 "-1" = (b+c)/(2*a+b+c)
# 3 "c" = (b+c)/2
# 4 "wb" = b+c
# 5 "r" = 2*b*c/((a+b+c)^2-2*b*c)
# 6 "I" = log(2*a+b+c) - 2*a*log(2)/(2*a+b+c) - ((a+b)*log(a+b) + (a+c)*log(a+c)) / (2*a+b+c)
# 7 "e" = exp(log(2*a+b+c) - 2*a*log(2)/(2*a+b+c) - ((a+b)*log(a+b) + (a+c)*log(a+c)) / (2*a+b+c))-1
# 8 "t" = (b+c)/(2*a+b+c)
# 9 "me" = (b+c)/(2*a+b+c)
# 10 "j" = a/(a+b+c)
# 11 "sor" = 2*a/(2*a+b+c)
# 12 "m" = (2*a+b+c)*(b+c)/(a+b+c)
# 13 "-2" = pmin(b,c)/(pmax(b,c)+a)
# 14 "co" = (a*c+a*b+2*b*c)/(2*(a+b)*(a+c))
# 15 "cc" = (b+c)/(a+b+c)
# 16 "g" = (b+c)/(a+b+c)
# 17 "-3" = pmin(b,c)/(a+b+c)
# 18 "l" = (b+c)/2
# 19 "19" = 2*(b*c+1)/(a+b+c)/(a+b+c-1)
# 20 "hk" = (b+c)/(2*a+b+c)
# 21 "rlb" = a/(a+c)
# 22 "sim" = pmin(b,c)/(pmin(b,c)+a)
# 23 "gl" = 2*abs(b-c)/(2*a+b+c)
# 24 "z" = (log(2)-log(2*a+b+c)+log(a+b+c))/log(2)

# calculating dissimilarity
j.dissimilarity <- betadiver(community.matrix, "g")

# 2 habitats
mod.habitat <- with(enviro.data, betadisper(j.dissimilarity, Habitat))
TukeyHSD.habitat <- TukeyHSD(mod.habitat)$group

# 4 climtes
mod.climate <- with(enviro.data, betadisper(j.dissimilarity, Climate))
TukeyHSD.climate <- TukeyHSD(mod.climate)$group

# 8 vegetations
mod.vegetation <- with(enviro.data, betadisper(j.dissimilarity, Vegetation))
TukeyHSD.vegetation <- TukeyHSD(mod.vegetation)$group

# save results
TukeyHSD.all <- rbind(TukeyHSD.habitat, TukeyHSD.climate, TukeyHSD.vegetation)
write.csv(TukeyHSD.all, file = "TukeyTest_all_17Jan2018.csv")

# combine all plots
layout(matrix(c(1,2,3,3),2,2, byrow=T))
boxplot(mod.habitat, ann=T, axes=T, col=c("grey30", "black"), border=par("fg"),
        las=1, pars=list(boxwex=0.2), cex.axis=1, ylim=c(0,1))
abline(h=0.6, lty=3)
boxplot(mod.climate,  ann=T, axes=T, col=c("grey20","dodgerblue3","green3","darkorange"), border=par("fg"),
        las=1, pars=list(boxwex=0.4), cex.axis=1, ylim=c(0,1))
abline(h=0.6, lty=3)
boxplot(mod.vegetation, ann=T, axes=T, col=c("darkviolet", "khaki", "goldenrod", "firebrick1", "darkgoldenrod1", "green3", "red", "orangered"),
        border=par("fg"),
        las=1, pars=list(boxwex=0.3), cex.axis=1, ylim=c(0,1))
abline(h=0.6, lty=3)


