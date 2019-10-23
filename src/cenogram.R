
# Ecological diversity analyses codes originally written by Meng Chen during the doctoral dissertation.

# 10July2017
# reorganized on 23Oct2019

######################################################
# Cenogram of body size of mammals of each community #
######################################################

library(vegan)
library(MASS)
library(permute)
library(calibrate)

# set up the directory
setwd("your directory")

### Using the average body size value
masterData <- read.csv(file="body size data",sep=",", header=T)

### number of order in the dataset
number.order <- length(unique(masterData$Order))
### number of family in the dataset
number.family <- length(unique(masterData$Family))
### number of species in the dataset
number.species <- length(unique(masterData$Species))

## data of different environments
bodysize.open.data <- masterData[masterData$Habitat == "Open",]
bodysize.close.data <- masterData[masterData$Habitat == "Close",]
bodysize.tropical.data <- masterData[masterData$Climate == "Tropical",]
bodysize.arid.data <- masterData[masterData$Climate == "Arid",]
bodysize.temperate.data <- masterData[masterData$Climate == "Temperate",]
bodysize.cold.data <- masterData[masterData$Climate == "Cold",]
bodysize.tropical.forest.data <- masterData[masterData$Vegetation == "Tropical rainforest",]
bodysize.tropical.seasonal.forest.data <- masterData[masterData$Vegetation == "Tropical seasonal forest",]
bodysize.savanna.data <- masterData[masterData$Vegetation == "Savanna",]
bodysize.grassland.data <- masterData[masterData$Vegetation == "Grassland",]
bodysize.shrubland.data <- masterData[masterData$Vegetation == "Shrubland",]
bodysize.desert.data <- masterData[masterData$Vegetation == "Desert",]
bodysize.temperate.forest.data <- masterData[masterData$Vegetation == "Temperate forest",]
bodysize.boreal.forest.data <- masterData[masterData$Vegetation == "Boreal forest",]

## all data as array
bodysize.data.array <- list(bodysize.open.data[,9], ## place holder
                            bodysize.open.data[,9],
                            bodysize.close.data[,9],
                            bodysize.close.data[,9],  ## place holder
                            bodysize.tropical.data[,9],
                            bodysize.arid.data[,9],
                            bodysize.temperate.data[,9],
                            bodysize.cold.data[,9],
                            bodysize.tropical.forest.data[,9],
                            bodysize.tropical.seasonal.forest.data[,9],
                            bodysize.savanna.data[,9],
                            bodysize.grassland.data[,9],
                            bodysize.shrubland.data[,9],
                            bodysize.desert.data[,9],
                            bodysize.temperate.forest.data[,9],
                            bodysize.boreal.forest.data[,9])
# number of the cenograms
number.cenogram <- length(bodysize.data.array)
# color array
cols.array <- c("black", "black", "grey30", "grey30",
                "darkorange", "grey20", "green3", "dodgerblue3",
                "red", "orangered", "firebrick1",
                "goldenrod", "darkgoldenrod1", "khaki", "green3", "darkviolet")
# enviromental names
diff.enviro <- c("","Open", "Close", "",
                 "Tropical", "Arid", "Temperate", "Cold", 
                 "Tropical rainforest", "Tropical seasonal forest", "Savanna", 
                 "Grassland", "Shrubland", "Desert", "Temperate forest", "Boreal forest")

par(mfcol=c(4,4),
    oma=c(3,3,1,1), 
    mar=c(2,2,1,1), 
    mgp=c(1,0.6,0),
    tck=-0.02)

for (i in 1:number.cenogram) {
  bs.sort <- sort(bodysize.data.array[[i]], decreasing=T)
  x = 1:length(bs.sort)
  y = log2(bs.sort)
  plot(x=x, y=y, ylim=c(0,14),
       las=1,
       pch=22,
       cex=1,
       bty="n",
       ann="F", 
       col=cols.array[i])
  mtext(outer=F,side=3, line=-1, at=length(bs.sort)/2, text=diff.enviro[i],cex=0.8)
}

