
# Resampling analyses codes originally written by Meng Chen during the doctoral dissertation.

# reorganized on 23Oct2019

################################################################
# The analyses here are focusing on the bootstrap resampling   #
# about communities in different environments in order to test #
# whether the sample size has a big impact on the results      #
# Author: Meng Chen                        Date: 12Jan2015     #
################################################################

# libraries used in the data analyses
library(permute)
library(lattice)
library(vegan)
library(MASS)
library(boot)
library(sp)
library(tidyr)
# set up the directory
setwd("your directory")

#################################################
# Parameters and Functions used in the analyses #
#################################################

# Parameters
habitat.pair <- c("open vs close")
climate.pair <- c("tropical vs arid", "tropical vs temperate", "tropical vs cold",
                  "arid vs temperate", "arid vs cold", "temperate vs cold")
vegetation.pair <- c("tropical rain forest vs tropical seasonal forest",
                     "tropical rain forest vs savanna",
                     "tropical rain forest vs grassland",
                     "tropical rain forest vs shrubland",
                     "tropical rain forest vs desert",
                     "tropical rain forest vs temperate forest",
                     "trapical rain forest vs boreal forest",
                     "tropical seasoanl forest vs savanna",
                     "tropical seasoanl forest vs grassland",
                     "tropical seasonal forest vs shrubland",
                     "tropical seasonal forest vs desert",
                     "tropical seasonal forest vs temperate forest",
                     "tropical seasonal forest vs boreal forest",
                     "savanna vs grassland",
                     "savanna vs shrubland",
                     "savanna vs desert",
                     "savanna vs temperate forest",
                     "savanna vs boreal forest",
                     "grassland vs shrubland",
                     "grassland vs desert",
                     "grassalnd vs temperate forest",
                     "grassland vs boreal forest",
                     "shrubland vs desert",
                     "shrubland vs temperate forest",
                     "shrubland vs boreal forest",
                     "desert vs temperate forest",
                     "desert vs boreal forest",
                     "temperate forest vs boreal forest")

# Resample function for speceis ecotypes most repeatedly
resample.func <- function(input.data, n.species) {
  set.seed(1234)
  resample.data <- sample(input.data[,1], 1000, replace=T)
  table.data <- table(resample.data)
  order.data <- table.data[order(table.data, decreasing=T)][1:n.species]
  ## ordered data becaome data frame with two columns
  data.frame.resample <- data.frame(order.data)
  ## separate the BS, Diet, and Locomotor
  resample.data.type <- cbind(read.fwf(file=textConnection(as.character(data.frame.resample[,1])),
                                       widths=c(1,1,1), colClasses="character",
                                       col.names=c("BS.Rank","Diet.Rank","Locomotor.Rank")),
                              data.frame.resample[-1])
   return(resample.data.type)
}

# function for resampled ERich based on average No of species
resampled.ERich.func <- function(input.data) {
  set.seed(1234)
  resampled.ERich <- data.frame()
  resampled.data <- sample(input.data[,1], 6, replace=T)
  temp.resampled.data <- unique(resampled.data)
  resampled.ERich <- length(temp.resampled.data)
  resampled.ERich
}

# function to calculate mean of resampled number of speceis in each evironmental type
mean.species.func <- function(input.data) {
  set.seed(1234)
  resampled.no.species <- sample(table(input.data[,1]), 1000, replace=T)
  mean.result <- mean(resampled.no.species)
  return(mean.result)
}

# writing table function
write.func <- function(data, filename) {
    write.table(data[,-4],
                file = filename,
                sep=",",
                row.names=F,
                col.names=F)
}

# Bootstrapping the specific number of species
bs.species.func <- function (input.data, n.species) {
  start <- Sys.time()
  mean.eds <- data.frame()
  set.seed(1234)
  for (t in 1:1000)
  {
    eds.community <- data.frame()
    data.temp <- input.data[sample(nrow(input.data), n.species, replace=T),]
    for (u in 1:(n.species-1))
    {
      for (v in u:(n.species-1))
      {
        results <- abs(data.temp[u,]-data.temp[v+1,])
        eds.temp <- rowSums(results)
        eds.community <- rbind(eds.community,c(eds.temp))
      }
    }
    mean.eds <- rbind(mean.eds,c(mean(eds.community[,1])))
  }
  print(Sys.time()-start)
  return(mean.eds)
}

# Histgram plot
hist.func <- function(data, bks, color) {
    h <- hist(data,
              breaks=bks,
              xlim=c(0,8),
              ylim=c(0,120),
              ann=F,
              col=color)
    xfit <- seq(min(data), max(data), length=40)
    yfit <- dnorm(xfit,mean=mean(data),sd=sd(data))
    yfit <- yfit*diff(h$mids[1:2])*length(data)
    lines(xfit, yfit, col="blue", lwd=1)
    abline(v=quantile(data, prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
}

#################################################################
# Bootstrapping for ecological function combination to identify #
# signture ecospace occupatons in each environmental type       #
#################################################################

# original extant small-bodied mammaian communities
mcm.habitat <- read.csv(file="EcoHabitatAll_6June2017.csv",
                        sep=",",
                        header=T)
## identify the russia park, four chinese sites, wetlands, and overgrassed
rownames(mcm.habitat[(mcm.habitat$Community.No == 68),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 72),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 102),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 105),])

# remove the unwanted communities
mcm.habitat <- mcm.habitat[-c(546:573, 821:823, 830:831),]

# Generate the ecotype of each species
mcm.habitat$Eco.type <- with(mcm.habitat, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

# 2 habitats
mcm.close <- mcm.habitat[mcm.habitat$Habitat == "Close",]
mcm.open <- mcm.habitat[mcm.habitat$Habitat == "Open",]

# 4 different climates
mcm.tropical.climate <- mcm.habitat[mcm.habitat$Climate == "Tropical",]
mcm.arid.climate <- mcm.habitat[mcm.habitat$Climate == "Arid",]
mcm.temperate.climate <- mcm.habitat[mcm.habitat$Climate== "Temperate",]
mcm.cold.climate <- mcm.habitat[mcm.habitat$Climate == "Cold",]

# 8 vegetation types
mcm.tropical.rainforest <- mcm.habitat[mcm.habitat$Vegetation == "Tropical rainforest",]
mcm.tropical.seasonal.forest <- mcm.habitat[mcm.habitat$Vegetation == "Tropical seasonal forest",]
mcm.savanna <- mcm.habitat[mcm.habitat$Vegetation == "Savanna",]
mcm.grassland <- mcm.habitat[mcm.habitat$Vegetation == "Grassland",]
mcm.shrubland <- mcm.habitat[mcm.habitat$Vegetation == "Shrubland",]
mcm.desert <- mcm.habitat[mcm.habitat$Vegetation == "Desert",]
mcm.temperate.forest <- mcm.habitat[mcm.habitat$Vegetation == "Temperate forest",]
mcm.boreal.forest <- mcm.habitat[mcm.habitat$Vegetation== "Boreal forest",]

# indicating how many average speceise of a community in each environment
n.species.close <- mean.species.func(mcm.close)
n.species.open <- mean.species.func(mcm.open)
n.species.tropical.climate <- mean.species.func(mcm.tropical.climate)
n.species.arid.climate  <- mean.species.func(mcm.arid.climate)
n.species.temperate.climate  <- mean.species.func(mcm.temperate.climate)
n.species.cold.climate  <- mean.species.func(mcm.cold.climate)
n.species.tropical.rainforest <- mean.species.func(mcm.tropical.rainforest)
n.species.tropical.seasonal.forest <- mean.species.func(mcm.tropical.seasonal.forest)
n.species.savanna <- mean.species.func(mcm.savanna)
n.species.grassland <- mean.species.func(mcm.grassland)
n.species.shrubland <- mean.species.func(mcm.shrubland)
n.species.desert <- mean.species.func(mcm.desert)
n.species.temperate.forest <- mean.species.func(mcm.temperate.forest)
n.species.boreal.forest <- mean.species.func(mcm.boreal.forest)

# the data from the Ecospace occpuation results
cl.type <- data.frame(mcm.close$Eco.type) # close habitat
op.type <- data.frame(mcm.open$Eco.type)## open habitat
tc.type <- data.frame(mcm.tropical.climate$Eco.type)## tropical climate
ac.type <- data.frame(mcm.arid.climate$Eco.type) ## arid climate
ttc.type <- data.frame(mcm.temperate.climate$Eco.type) ## temperate climate
cc.type <- data.frame(mcm.cold.climate$Eco.type)## cold climate
trf.type <- data.frame(mcm.tropical.rainforest$Eco.type) ## tropical rianforest
tsf.type <- data.frame(mcm.tropical.seasonal.forest$Eco.type) ## tropical seasonal forest
sv.type <- data.frame(mcm.savanna$Eco.type)## savanna
gl.type <- data.frame(mcm.grassland$Eco.type)  ## grassland
sl.type <- data.frame(mcm.shrubland$Eco.type)## shrubland
dt.type <- data.frame(mcm.desert$Eco.type) ## desert
ttf.type <- data.frame(mcm.temperate.forest$Eco.type)## temperate forest
bf.type <- data.frame(mcm.boreal.forest$Eco.type)## boreal forest

# two habiats
cl.resample.type <- resample.func(cl.type, n.species.close)
write.func(cl.resample.type, "sig.close_16Jan2018.csv")
op.resample.type <- resample.func(op.type, n.species.open)
write.func(op.resample.type, "sig.open_16Jan2018.csv")

# four climates
tc.resample.type <- resample.func(tc.type, n.species.tropical.climate)
write.func(tc.resample.type, "sig.tropical.climate_16Jan2018.csv")
ac.resample.type <- resample.func(ac.type, n.species.arid.climate)
write.func(ac.resample.type, "sig.arid.climate_16Jan2018.csv")
ttc.resample.type <- resample.func(ttc.type, n.species.temperate.climate)
write.func(ttc.resample.type, "sig.temperate.climate_16Jan2018.csv")
cc.resample.type <- resample.func(cc.type, n.species.cold.climate)
write.func(cc.resample.type, "sig.cold.climate_16Jan2018.csv")

# eight vegetations
trf.resample.type <- resample.func(trf.type, n.species.tropical.rainforest)
write.func(trf.resample.type, "sig.tropical.rainforest_16Jan2018.csv")
tsf.resample.type <- resample.func(tsf.type, n.species.tropical.seasonal.forest)
write.func(tsf.resample.type, "sig.tropical.seasonal.forest_16Jan2018.csv")
sv.resample.type <- resample.func(sv.type, n.species.savanna)
write.func(sv.resample.type, "sig.savanna_16Jan2018.csv")
gl.resample.type <- resample.func(gl.type, n.species.grassland)
write.func(gl.resample.type, "sig.grassland_16Jan2018.csv")
sl.resample.type <- resample.func(sl.type, n.species.shrubland)
write.func(sl.resample.type, "sig.shrubland_16Jan2018.csv")
dt.resample.type <- resample.func(dt.type, n.species.desert)
write.func(dt.resample.type, "sig.desert_16Jan2018.csv", sep=",")
ttf.resample.type <- resample.func(ttf.type, n.species.temperate.forest)
write.func(ttf.resample.type, "sig.temperate.forest_16Jan2018.csv")
bf.resample.type <- resample.func(bf.type, n.species.boreal.forest)
write.func(bf.resample.type, "sig.boreal.forest_16Jan2018.csv")

# convert the resample dataset to a table
sig.mastersheet <- read.csv(file="Sig_MaterSheet_Resampled_16Jan2018.csv", sep=",")
sig.mastersheet <- sig.mastersheet[,c(1,5)]
write.table(table(sig.mastersheet), file="Sig_MaterSheet_Summary_Resampled_16Jan2018.csv", sep=",")

####################################################################
# Bootstrapping for designated numbers of species with a community #
####################################################################
# Caution:
# Don't run it again because the resampling data has been saved already.
# If you accidentally run the resamppling code, the plots and distribution of
# may be shifted a little bit. However, all plots are needed to be revised.

###########################################
# 5 Species in each simulated communities #
###########################################

# 2 habitats
# Closed habitat
cl.results <- bs.species.func(mcm.close[,3:5], 5)
# Open habitat
op.results <- bs.species.func(mcm.open[,3:5], 5)

# 4 different climates
# Tropical climate
tc.results <- bs.species.func(mcm.tropical.climate[,3:5], 5)
# Arid climate
ac.results <- bs.species.func(mcm.arid.climate[,3:5], 5)
# Temperate climate
ttc.results <- bs.species.func(mcm.temperate.climate[,3:5], 5)
# Cold climate
cc.results <- bs.species.func(mcm.cold.climate[,3:5], 5)

# 8 vegetations
# Tropical rainforest
trf.results <- bs.species.func(mcm.tropical.rainforest[,3:5], 5)
# Tropical seasonal forest
tsf.results <- bs.species.func(mcm.tropical.seasonal.forest[,3:5], 5)
# Savanna
sv.results <- bs.species.func(mcm.savanna[,3:5], 5)
# grassland
gl.results <- bs.species.func(mcm.grassland[,3:5], 5)
# Shrubland
sl.results <- bs.species.func(mcm.shrubland[,3:5], 5)
# Desert
dt.results <- bs.species.func(mcm.desert[,3:5], 5)
# Temperate forest
ttf.results <- bs.species.func(mcm.temperate.forest[,3:5], 5)
# Boreal forest
bf.results <- bs.species.func(mcm.boreal.forest[,3:5], 5)

# Mean and SD of all bootstrapped values
bootstrapped.results.all.5p <- matrix(c(cl.results[,1], op.results[,1],
                                        tc.results[,1], ac.results[,1],
                                        ttc.results[,1], cc.results[,1],
                                        trf.results[,1], tsf.results[,1],
                                        sv.results[,1], gl.results[,1],
                                        sl.results[,1], dt.results[,1],
                                        ttf.results[,1], bf.results[,1]),
                                      nrow = 1000,
                                      ncol = 14,
                                      byrow = F)
colnames(bootstrapped.results.all.5p) <- c("Close", "Open", "Tropical","Arid",
                                           "Temperate", "Cold", "Tropical rainforest",
                                           "Tropical seasonal forest", "Savanna",
                                           "Grassland", "Shrubland", "Desert",
                                           "Temperate forest", "Boreal forest")
mean.bootstrapped.5p <- apply(bootstrapped.results.all.5p, 2, mean)
sd.bootstrapped.5p <- apply(bootstrapped.results.all.5p, 2, sd)
qt.bootstrapped.5p <- apply(bootstrapped.results.all.5p, 2, quantile, prob=c(0.025, 0.975))

## save the data
bootstropped.results.stat.5p <- t(rbind(t(mean.bootstrapped.5p), t(sd.bootstrapped.5p)))
colnames(bootstropped.results.stat.5p) <- c("Mean", "sd")
bootstropped.results.stat.5p <- cbind(bootstropped.results.stat.5p, t(qt.bootstrapped.5p))
write.csv(bootstropped.results.stat.5p, file="Bootstrapped_statistic_5species_16Jan2018.csv")

####################
# Pair-wise t test #
####################

# two habitats
habitat.results.test.5p <- data.frame()
habitat.results.test.5p <- t.test(bootstrapped.results.all.5p[,1],
                                  bootstrapped.results.all.5p[,2])
habitat.results.test.5p <- rbind(c(habitat.results.test.5p$statistic,
                                   habitat.results.test.5p$parameter,
                                   habitat.results.test.5p$p.value))
colnames(habitat.results.test.5p) <- c("t","df","P value")
rownames(habitat.results.test.5p) <- habitat.pair

# four climates
climate.results.test.5p <- data.frame()
for (z in 3:5)
{
  for (h in z:5) {
    temp.result.test.5p <- t.test(bootstrapped.results.all.5p[,z],
                                  bootstrapped.results.all.5p[, h+1])
    climate.results.test.5p <- rbind(climate.results.test.5p,
                                     c(temp.result.test.5p$statistic,
                                       temp.result.test.5p$parameter,
                                       temp.result.test.5p$p.value))
  }
}
colnames(climate.results.test.5p)<-c("t","df","P value")
rownames(climate.results.test.5p) <- climate.pair

#### 8 vegetations
vegetation.results.test.5p<-data.frame()
for (w in 7:13)
{
  for (q in w:13) {
    temp.result.test.5p<-t.test(bootstrapped.results.all.5p[,w],
                                bootstrapped.results.all.5p[,q+1])
    vegetation.results.test.5p<-rbind(vegetation.results.test.5p,
                                      c(temp.result.test.5p$statistic,
                                        temp.result.test.5p$parameter,
                                        temp.result.test.5p$p.value))
  }
}
colnames(vegetation.results.test.5p)<-c("t","df","P value")
rownames(vegetation.results.test.5p) <- vegetation.pair

## combine all results
t.test.results.all.5species <- rbind(habitat.results.test.5p,
                                     climate.results.test.5p,
                                     vegetation.results.test.5p)
write.csv(t.test.results.all.5species,
          file = "Bootstrapped_statistic_table_5species_16Jan2018.csv")

## histgram plots
par(mfcol=c(8,3),
    oma=c(1,1,0,1),
    mar=c(2,2,1,0),
    mgp=c(1,0.6,0),
    tck=-0.02)

# 2 habitats
# filling the empty space of the composite plot
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
# closed habitat
hist.func(cl.results[,1], 40, "grey30")
# open habitat
hist.func(op.results[,1], 40, "black")
# filling the empty space of the composite plot
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)

# 4 climates
# filling the empty space of the composite plot
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
# tropical climte
hist.func(tc.results[,1], 30, "darkorange")
# arid climate
hist.func(ac.results[,1], 40, "grey20")
# temperate climate
hist.func(ttc.results[,1], 40, "green3")
# cold climate
hist.func(cc.results[,1], 40, "dodgerblue3")
# filling the empty space of the composite plot
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)

# 8 vegetations
# tropical rainforest
hist.func(trf.results[,1], 30, "red")
# tropical seasonal forest
hist.func(tsf.results[,1], 40, "orangered")
# savanna
hist.func(sv.results[,1], 40, "firebrick1")
# grassland
hist.func(gl.results[,1], 40, "goldenrod)
# shrubland
hist.func(sl.results[,1], 40, "darkgoldenrod1")
# desert
hist.func(dt.results[,1], 30, "khaki")
# temperate forest
hist.func(ttf.results[,1], 40, "green3")
# boreal forest
hist.func(bf.results[,1], 40, "darkviolet")

###########################################
# 6 Species in each simulated communities #
###########################################

# 2 habitats
# Closed habitat
cl.results <- bs.species.func(mcm.close[,3:5], 6)
# Open habitat
op.results <- bs.species.func(mcm.open[,3:5], 6)

# 4 different climates
# Tropical climate
tc.results <- bs.species.func(mcm.tropical.climate[,3:5], 6)
# Arid climate
ac.results <- bs.species.func(mcm.arid.climate[,3:5], 6)
# Temperate climate
ttc.results <- bs.species.func(mcm.temperate.climate[,3:5], 6)
## Cold climate
cc.results <- bs.species.func(mcm.cold.climate[,3:5], 6)

# 8 vegetations
# Tropical rainforest
trf.results <- bs.species.func(mcm.tropical.rainforest[,3:5], 6)
# Tropical seasonal forest
tsf.results <- bs.species.func(mcm.tropical.seasonal.forest[,3:5], 6)
# Savanna
sv.results <- bs.species.func(mcm.savanna[,3:5], 6)
# Grassland
gl.results <- bs.species.func(mcm.grassland[,3:5], 6)
# Shrubland
sl.results <- bs.species.func(mcm.shrubland[,3:5], 6)
# Desert
dt.results <- bs.species.func(mcm.desert[,3:5], 6)
# Temperate forest
ttf.results <- bs.species.func(mcm.temperate.forest[,3:5], 6)
# Boreal forest
bf.results <- bs.species.func(mcm.boreal.forest[,3:5], 6)

## Mean and SD of all bootstrapped values
bootstrapped.results.all.6p <- matrix(c(cl.results[,1], op.results[,1],
                                        tc.results[,1], ac.results[,1],
                                        ttc.results[,1], cc.results[,1],
                                        trf.results[,1], tsf.results[,1],
                                        sv.results[,1], gl.results[,1],
                                        sl.results[,1], dt.results[,1],
                                        ttf.results[,1], bf.results[,1]),
                                      nrow = 1000, ncol = 14, byrow = F)

colnames(bootstrapped.results.all.6p) <- c("Close", "Open", "Tropical","Arid",
                                           "Temperate", "Cold", "Tropical rainforest",
                                           "Tropical seasonal forest", "Savanna",
                                           "Grassland", "Shrubland", "Desert",
                                           "Temperate forest", "Boreal forest")
mean.bootstrapped.6p <- apply(bootstrapped.results.all.6p, 2, mean)
sd.bootstrapped.6p <- apply(bootstrapped.results.all.6p, 2, sd)
qt.bootstrapped.6p <- apply(bootstrapped.results.all.6p, 2, quantile, prob=c(0.025, 0.975))

## manage the data
bootstropped.results.stat.6p <- t(rbind(t(mean.bootstrapped.6p), t(sd.bootstrapped.6p)))
colnames(bootstropped.results.stat.6p) <- c("Mean", "sd")
bootstropped.results.stat.6p <- cbind(bootstropped.results.stat.6p, t(qt.bootstrapped.6p))
write.csv(bootstropped.results.stat.6p, file="Bootstrapped_statistic_6species_16Jan2018.csv")

####################
# Pair-wise t test #
####################

# two habitats
habitat.results.test.6p <- data.frame()
habitat.results.test.6p <- t.test(bootstrapped.results.all.6p[,1],
                                  bootstrapped.results.all.6p[,2])
habitat.results.test.6p <- rbind(c(habitat.results.test.6p$statistic,
                                   habitat.results.test.6p$parameter,
                                   habitat.results.test.6p$p.value))
colnames(habitat.results.test.6p) <- c("t","df","P value")
rownames(habitat.results.test.6p) <- habitat.pair

# four climates
climate.results.test.6p <- data.frame()
for (z in 3:5)
{
  for (h in z:5) {
    temp.result.test.6p <- t.test(bootstrapped.results.all.6p[,z],
                                  bootstrapped.results.all.6p[, h+1])
    climate.results.test.6p <- rbind(climate.results.test.6p,
                                     c(temp.result.test.6p$statistic,
                                       temp.result.test.6p$parameter,
                                       temp.result.test.6p$p.value))
  }
}
colnames(climate.results.test.6p) <- c("t","df","P value")
rownames(climate.results.test.6p) <- climate.pair

# eight vegetations
vegetation.results.test.6p <- data.frame()
for (w in 7:13)
{
  for (q in w:13) {
    temp.result.test.6p <- t.test(bootstrapped.results.all.6p[,w],
                                  bootstrapped.results.all.6p[, q+1])
    vegetation.results.test.6p <- rbind(vegetation.results.test.6p,
                                        c(temp.result.test.6p$statistic,
                                          temp.result.test.6p$parameter,
                                          temp.result.test.6p$p.value))
  }
}
colnames(vegetation.results.test.6p) <- c("t","df","P value")
rownames(vegetation.results.test.6p) <- vegetation.pair

## combine all results
t.test.results.all.6species <- rbind(habitat.results.test.6p,
                                     climate.results.test.6p,
                                     vegetation.results.test.6p)
write.csv(t.test.results.all.6species, file = "Bootstrapped_statistic_table_6species_16Jan2018.csv")

## histgram plots
par(mfcol=c(8,3),
    oma=c(1,1,0,1),
    mar=c(2,2,1,0),
    mgp=c(1,0.6,0),
    tck=-0.02)

# 2 habitats
# filling the empty space of the composite plot
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
# closed habitat
hist.func(cl.results[,1], 40, "grey30")
# open habitat
hist.func(op.results[,1], 30, "black")
# filling the empty space of the composite plot
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)

# 4 climates
# filling the empty space of the composite plot
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
# tropical climte
hist.func(tc.results[,1], 30, "darkorange")
# arid climate
hist.func(ac.results[,1], 40, "grey20")
# temperate climate
hist.func(ttc.results[,1], 40, "green3")
# cold climate
hist.func(cc.results[,1], 30, "dodgerblue3")
# filling the empty space of the composite plot
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)

# 8 vegetations
# tropical rainforest
hist.func(trf.results[,1], 30, "red")
# tropical seasonal forest
hist.func(tsf.results[,1], 40, "orangered")
# savanna
hist.func(sv.results[,1], 40, "firebrick1")
# grassland
hist.func(gl.results[,1], 30, "goldenrod)
# shrubland
hist.func(sl.results[,1], 30, "darkgoldenrod1")
# desert
hist.func(dt.results[,1], 30, "khaki")
# temperate forest
hist.func(ttf.results[,1], 30, "green3")
# boreal forest
hist.func(bf.results[,1], 30, "darkviolet")

############################################
# 10 Species in each simulated communities #
############################################

# 2 habitats
# Closed habitat
cl.results <- bs.species.func(mcm.close[,3:5], 10)
# Open
op.results <- bs.species.func(mcm.open[,3:5], 10)

# 4 different climates
# Tropical climate
tc.results <- bs.species.func(mcm.tropical.climate[,3:5], 10)
# Arid climate
ac.results <- bs.species.func(mcm.arid.climate[,3:5], 10)
# Temperate climate
ttc.results <- bs.species.func(mcm.temperate.climate[,3:5], 10)
## Cold climate
cc.results <- bs.species.func(mcm.cold.climate[,3:5], 10)

# 8 vegetations
# Tropical rainforest
trf.results <- bs.species.func(mcm.tropical.rainforest[,3:5], 10)
# Tropical seasonal forest
tsf.results <- bs.species.func(mcm.tropical.seasonal.forest[,3:5], 10)
# Savanna
sv.results <- bs.species.func(mcm.savanna[,3:5], 10)
# Grassland
gl.results <- bs.species.func(mcm.grassland[,3:5], 10)
# Shrubland
sl.results <- bs.species.func(mcm.shrubland[,3:5], 10)
# Desert
dt.results <- bs.species.func(mcm.desert[,3:5], 10)
# Temperate forest
ttf.results <- bs.species.func(mcm.temperate.forest[,3:5], 10)
# Boreal forest
bf.results <- bs.species.func(mcm.boreal.forest[,3:5], 10)

## Mean and SD of all bootstrapped values
bootstrapped.results.all.10p <- matrix(c(cl.results[,1], op.results[,1],
                                         tc.results[,1], ac.results[,1],
                                         ttc.results[,1], cc.results[,1],
                                         trf.results[,1], tsf.results[,1],
                                         sv.results[,1], gl.results[,1],
                                         sl.results[,1], dt.results[,1],
                                         ttf.results[,1], bf.results[,1]),
                                       nrow = 1000, ncol = 14, byrow = F)

colnames(bootstrapped.results.all.10p) <- c("Close", "Open", "Tropical","Arid",
                                            "Temperate", "Cold", "Tropical rainforest",
                                            "Tropical seasonal forest", "Savanna",
                                            "Grassland", "Shrubland", "Desert",
                                            "Temperate forest", "Boreal forest")
mean.bootstrapped.10p <- apply(bootstrapped.results.all.10p, 2, mean)
sd.bootstrapped.10p <- apply(bootstrapped.results.all.10p, 2, sd)
qt.bootstrapped.10p <- apply(bootstrapped.results.all.10p, 2, quantile, prob=c(0.025, 0.975))

## manage the data
bootstropped.results.stat.10p <- t(rbind(t(mean.bootstrapped.10p), t(sd.bootstrapped.10p)))
colnames(bootstropped.results.stat.10p) <- c("Mean", "sd")
bootstropped.results.stat.10p <- cbind(bootstropped.results.stat.10p, t(qt.bootstrapped.10p))
write.csv(bootstropped.results.stat.10p, file="Bootstrapped_statistic_10species_16Jan2018.csv")

#####################
# Pair-wise t test  #
#####################

# 2 habitats
habitat.results.test.10p <- data.frame()
habitat.results.test.10p <- t.test(bootstrapped.results.all.10p[,1],
                                   bootstrapped.results.all.10p[,2])
habitat.results.test.10p <- rbind(c(habitat.results.test.10p$statistic,
                                    habitat.results.test.10p$parameter,
                                    habitat.results.test.10p$p.value))
colnames(habitat.results.test.10p) <- c("t","df","P value")
rownames(habitat.results.test.10p) <- habitat.pair

# 4 climate
climate.results.test.10p <- data.frame()
for (z in 3:5)
{
  for (h in z:5) {
    temp.result.test.10p <- t.test(bootstrapped.results.all.10p[,z],
                                   bootstrapped.results.all.10p[, h+1])
    climate.results.test.10p <- rbind(climate.results.test.10p,
                                      c(temp.result.test.10p$statistic,
                                        temp.result.test.10p$parameter,
                                        temp.result.test.10p$p.value))
  }
}
colnames(climate.results.test.10p) <- c("t","df","P value")
rownames(climate.results.test.10p) <- climate.pair

# 8 vegetations
vegetation.results.test.10p <- data.frame()
for (w in 7:13)
{
  for (q in w:13) {
    temp.result.test.10p <- t.test(bootstrapped.results.all.10p[,w],
                                   bootstrapped.results.all.10p[, q+1])
    vegetation.results.test.10p <- rbind(vegetation.results.test.10p,
                                         c(temp.result.test.10p$statistic,
                                           temp.result.test.10p$parameter,
                                           temp.result.test.10p$p.value))
  }
}
colnames(vegetation.results.test.10p) <- c("t","df","P value")
rownames(vegetation.results.test.10p) <- vegetation.pair

## combine all results
t.test.results.all.10species <- rbind(habitat.results.test.10p,
                                      climate.results.test.10p,
                                      vegetation.results.test.10p)
write.csv(t.test.results.all.10species, file = "Bootstrapped_statistic_table_10species_16Jan2018.csv")

## histgram plots
par(mfcol=c(8,3),
    oma=c(1,1,0,1),
    mar=c(2,2,1,0),
    mgp=c(1,0.6,0),
    tck=-0.02)

# 2 habitats
# filling the empty space of the composite plot
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
# closed habitat
hist.func(cl.results[,1], 20, "grey30")
# open habitat
hist.func(op.results[,1], 20, "back")
# filling the empty space of the composite plot
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F)

# 4 climates
# filling the empty space of the composite plot
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
# tropical climte
hist.func(tc.results[,1], 25, "darkorange")
# arid climate
hist.func(ac.results[,1], 30, "grey20")
# temperate climate
hist.func(ttc.results[,1], 20, "green3")
# cold climate
hist.func(cc.results[,1], 20, "dodgerblue3")
# filling the empty space of the composite plot
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F)

# 8 vegetations
# tropical rainforest
hist.func(trf.results[,1], 20, "red")
# tropical seasonal forest
hist.func(tsf.results[,1], 30, "orangered")
# savanna
hist.func(sv.results[,1], 30, "firebrick1")
# grassland
hist.func(gl.results[,1], 20, "goldenrod)
# shrubland
hist.func(sl.results[,1], 20, "darkgoldenrod1")
# desert
hist.func(dt.results[,1], 15, "khaki")
# temperate forest
hist.func(ttf.results[,1], 30, "green3")
# boreal forest
hist.func(bf.results[,1], 30, "darkviolet")
