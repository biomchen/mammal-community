
# ggplot2 codes originally written by Meng Chen for his doctoral dissertation.

# 12Jan2015
# reorganized on 23Oct2019

# libraries used in the data analyses
library(permute)
library(vegan)
library(MASS)
library(ggplot2)
library(gridExtra)

# set up your directory
setwd("your directory")

# Eco dataset
EDisp_ERich <- read.csv(file = "EDisp and ERich data", sep = ",", header = T)
DF1_TwoHabitats <- read.csv(file = "DF1s of 2 habitat data", sep = ",", header = T)
# combine the datasets
combined.eco.DF1 <- cbind(DF1_TwoHabitats, EDisp_ERich)
combined.eco.DF1 <- combined.eco.DF1[,-c(1, 5)]
# the number of species dataset
mcm.habitat<-read.csv(file="mammal community eco data", sep=",", header=T)
n.species.each.community <- table(mcm.habitat$Community.No)
n.species.each.community <- as.data.frame(n.species.each.community)
colnames(n.species.each.community) <- c("Community", "No_species")

# GPS coordinate data
# get data
coords.eco <- read.csv(file="GPS coordinate data", header=TRUE, sep = ",")
final.combined.eco.data <- cbind(coords.eco, combined.eco.DF1[,3:5])
# numbers of the communities in each habitat
n.close <- dim(coords.eco[coords.eco$Habitat=="Close",])[1]
n.open <- dim(coords.eco[coords.eco$Habitat=="Open",])[1]

# save the data for the manually add the GPS points
write.csv(final.combined.eco.data, file = "your file name")

###########################
# ERich vs No. of species #
###########################

# regresson
fit.ERich.NoSpecies <- lm(n.species.each.community$No_species~combined.eco.DF1$ERich)
# combined data
NoSpecies.ERich.data <- cbind(n.species.each.community$No_species, combined.eco.DF1$ERich)
colnames(NoSpecies.ERich.data) <- c("No_species", "ERich")
NoSpecies.ERich.data <- as.data.frame(NoSpecies.ERich.data)
# linear model plot
ggplot(data=NoSpecies.ERich.data, aes(x=No_species, y=ERich))+ geom_point(shape=1) + geom_smooth(method=lm)

#########################################################
# ploting EDisp and ERich against latitudinal gradients #
#########################################################

# get the data
Eco.data.latitude <- read.csv(file="latitudinal gradient data", header=T, sep =",")
# converted the data
Eco.data.latitude$abs.Latitude <- abs(Eco.data.latitude$Latitude)
Eco.data.latitude <- as.data.frame(Eco.data.latitude)

# all communities
# colors
cols.habitat=c(rep("grey30", n.close),  # close
               rep("black", n.open),  # open
               rep("grey30", 5))           
# symbols
sybs.habitat=c(rep(19, n.close),       # close   
               rep(1, n.open),         # open     
               rep(16,5))
# initial visualization for all communities
plot(DF1~Latitude, 
     pch = c(6, 19, 1)[as.numeric(Habitat)],  ## 6 for NA, 4 for close, 15 for open
     col= c("grey30", "grey20","black")[as.numeric(Habitat)],
     data = Eco.data.latitude,
     xlim = c(-50, 100),
     cex=2)
plot(Mean.EDisp~Latitude, 
     type = "h",
     data = Eco.data.latitude,
     xlim = c(-50, 100),
     ylim = c(0, 8),
     cex=1.5)
plot(ERich~Latitude, 
     type = "h",
     data = Eco.data.latitude,
     xlim = c(-50, 100),
     ylim = c(0, 30),
     cex=1.5)

# visualization using different enviromental categories
# two habitats
# EDisp
p1 <- ggplot(Eco.data.latitude, aes(Latitude, Mean.EDisp)) +
  geom_bar(stat="identity",aes(color=Habitat), size=0.5, alpha=0.5)  +
  scale_color_manual(values=c("grey30", "black")) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")
# ERich
p2 <- ggplot(Eco.data.latitude, aes(Latitude, ERich)) +
  geom_bar(stat="identity",aes(color=Habitat), size=0.5, alpha=0.5)  +
  scale_color_manual(values=c("grey30", "black")) + 
  coord_flip() +
  theme(legend.position="none")

# four climates
# EDisp
p3 <- ggplot(Eco.data.latitude, aes(Latitude, Mean.EDisp)) +
  geom_bar(stat="identity",aes(color=Climate), size=0.5, alpha=0.5) + 
  scale_color_manual(values=c("grey20", "dodgerblue3", "green3", "darkorange")) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")
# ERich
p4 <- ggplot(Eco.data.latitude, aes(Latitude, ERich)) +
  geom_bar(stat="identity",aes(color=Climate), size=0.5, alpha=0.5) + 
  scale_color_manual(values=c("grey20", "dodgerblue3", "green3", "darkorange")) +
  coord_flip() +
  theme(legend.position="none")

# eight vegetations
# EDisp
p5 <- ggplot(Eco.data.latitude, aes(Latitude, Mean.EDisp)) +
  geom_bar(stat="identity",aes(color=Vegetation), size=0.5, alpha=0.5) + 
  scale_color_manual(values=c("darkviolet", "khaki", "goldenrod", "firebrick1", "darkgoldenrod1", "green3", "red", "orangered")) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")
# ERich
p6 <- ggplot(Eco.data.latitude, aes(Latitude, ERich)) +
  geom_bar(stat="identity",aes(color=Vegetation), size=0.5, alpha=0.5) + 
  scale_color_manual(values=c("darkviolet", "khaki", "goldenrod", "firebrick1", "darkgoldenrod1", "green3", "red", "orangered")) +
  coord_flip() + 
  theme(legend.position="none")

# Mean Annual Temperature (MAT)
# EDisp
p7 <- ggplot(Eco.data.latitude, aes(Latitude, Mean.EDisp)) +
  geom_bar(stat="identity",aes(color=MAT), size=0.5, alpha=0.5) + 
  scale_colour_gradient(low="#4169E1" ,high="#EE5C42", breaks=c(24, 16, 8, 0)) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")
# ERich
p8 <- ggplot(Eco.data.latitude, aes(Latitude, ERich)) +
  geom_bar(stat="identity",aes(color=MAT), size=0.5, alpha=0.5) + 
  scale_colour_gradient(low="#4169E1" ,high="#EE5C42", breaks=c(24, 16, 8, 0)) +
  coord_flip()  +
  theme(legend.position="none")

# Mean Annual Precipitation (MAP)
# EDisp
p9 <- ggplot(Eco.data.latitude, aes(Latitude, Mean.EDisp)) +
  geom_bar(stat="identity",aes(color=MAP), size=0.5, alpha=0.2) + 
  scale_colour_gradient(low="#E69F00" ,high="#00FF7F", breaks=c(6000, 4800, 3600, 2400, 1200)) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")
# ERich
p10 <- ggplot(Eco.data.latitude, aes(Latitude, ERich)) +
  geom_bar(stat="identity",aes(color=MAP), size=0.5, alpha=0.2) + 
  scale_colour_gradient(low="#E69F00" ,high="#00FF7F", breaks=c(6000, 4800, 3600, 2400, 1200)) +
  coord_flip() +
  theme(legend.position="none")

# Elevation
# EDisp
p11 <- ggplot(Eco.data.latitude, aes(Latitude, Mean.EDisp)) +
  geom_bar(stat="identity",aes(color=Elevation), size=0.5, alpha=0.2) + 
  scale_color_gradient(low="grey80", high="grey10", breaks=c(3200, 2400, 1600, 800)) +
  scale_y_reverse() + 
  coord_flip() + 
  theme(legend.position="none")
# ERich
p12 <- ggplot(Eco.data.latitude, aes(Latitude, ERich)) +
  geom_bar(stat="identity",aes(color=Elevation), size=0.5, alpha=0.2) + 
  scale_color_gradient(low="grey80", high="grey10", breaks=c(3200, 2400, 1600, 800)) +
  coord_flip() +
  theme(legend.position="none")

# all plots together
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol=6, nrow=2)

####################################
# Latidutional gradients along DF1 #
####################################

# two habitats
p13 <- ggplot(Eco.data.latitude, aes(Latitude, DF1)) +
  geom_point(aes(color=Habitat, shape=Habitat), size=3, alpha=1) +
  scale_color_manual(values=c("red", "grey30", "black")) +  ## red is for NA
  scale_shape_manual(values=c(9, 19, 1)) +   ## 9 is for NA
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")

# four climates
p14 <- ggplot(Eco.data.latitude, aes(Latitude, DF1)) +
  geom_point(aes(color=Climate, shape=Climate), size=3, alpha=1) +
  scale_color_manual(values=c("red", "grey20","dodgerblue3","green3","darkorange")) +  ## red is for NA
  scale_shape_manual(values=c(9, 18, 15, 17, 16)) + ## 9 is for NA
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")

# eight vegetations
p15 <- ggplot(Eco.data.latitude, aes(Latitude, DF1)) +
  geom_point(aes(color=Vegetation, shape=Vegetation), size=3, alpha=1) +
  scale_color_manual(values=c("red", "darkviolet", "khaki", "goldenrod", "firebrick1", "darkgoldenrod1", "green3","red", "orangered")) +  ## red is for NA
  scale_shape_manual(values=c(9, 1, 18, 15, 2, 4, 0, 16, 11)) +  ## 9 is for NA
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")

# MAT
p16 <- ggplot(Eco.data.latitude, aes(Latitude, DF1)) +
  geom_point(aes(color=MAT), size=3, alpha=1) +
  scale_color_gradient(low="#4169E1" ,high="#EE5C42", breaks=c(24,16,8, 0)) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")

# MAP
p17 <- ggplot(Eco.data.latitude, aes(Latitude, DF1)) +
  geom_point(aes(color=MAP), size=3, alpha=1) +
  scale_colour_gradient(low="#E69F00" ,high="#00FF7F", breaks=c(6000,4800,3600,2400,1200)) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")

# Elevation
p18 <- ggplot(Eco.data.latitude, aes(Latitude, DF1)) +
  geom_point(aes(color=Elevation), size=3, alpha=1) +
  scale_color_gradient(low="grey80", high="grey10", breaks=c(3200, 2400, 1600, 800)) +
  scale_y_reverse() + 
  coord_flip() +
  theme(legend.position="none")

# all plots together
grid.arrange(p13, p14, p15, p16, p17, p18, ncol=3, nrow=2)
