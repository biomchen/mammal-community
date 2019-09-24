
### Ecological diversity analyses codes originally written by Meng Chen during the doctoral dissertation
### Reorder method is using climate region based on climate codes

### The purpose of the analyses is to aim to calculate how many different combination of three ecological
### paramters in each mammal community to decipher the ecological differences among them

#########################################################################################################
# The analyses here are focusing on the bootstrap resampling about the communities in different climate #
# environment in order to test whether the sample size has a big impact on the results                  #
# Meng Chen                                                           Date: 12Jan2016                   #
#########################################################################################################

# Possible libraries used in the data analyses
library(permute)
library(lattice)
library(vegan)
library(MASS)
library(boot)
library(sp)
library(parallel)
library(tidyr)

setwd("/Users/mengchen/Documents/GitHub/2018_Mamm_Comm_Evol_Ecol⁩/data⁩")

#########################################################################################################################
# Bootstrapping for ecological function combination to identify signture ecospace occupatons in each environmental type #
#########################################################################################################################

# original extant small-bodied mammaian communities
mcm.habitat<-read.csv(file="EcoHabitatAll_6June2017.csv", sep=",", header=T)
## identify the russia park, four chinese sites, wetlands, and overgrassed
rownames(mcm.habitat[(mcm.habitat$Community.No == 68),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 72),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 102),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 105),])

# remove the unwanted communities
mcm.habitat<- mcm.habitat[-c(546:573, 821:823, 830:831),]

# Generate the ecotype of each species
mcm.habitat$Eco.type <- with(mcm.habitat, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

# 2 habitats
mcm.close<-mcm.habitat[mcm.habitat$Habitat == "Close",]
mcm.open<-mcm.habitat[mcm.habitat$Habitat == "Open",]

# 4 different climates
mcm.tropical.climate<-mcm.habitat[mcm.habitat$Climate == "Tropical",]
mcm.arid.climate<-mcm.habitat[mcm.habitat$Climate == "Arid",]
mcm.temperate.climate<-mcm.habitat[mcm.habitat$Climate== "Temperate",]
mcm.cold.climate<-mcm.habitat[mcm.habitat$Climate == "Cold",]

# 8 vegetation types
mcm.tropical.rainforest<-mcm.habitat[mcm.habitat$Vegetation == "Tropical rainforest",]
mcm.tropical.seasonal.forest<-mcm.habitat[mcm.habitat$Vegetation == "Tropical seasonal forest",]
mcm.savanna<-mcm.habitat[mcm.habitat$Vegetation == "Savanna",]
mcm.grassland<-mcm.habitat[mcm.habitat$Vegetation == "Grassland",]
mcm.shrubland<-mcm.habitat[mcm.habitat$Vegetation == "Shrubland",]
mcm.desert<-mcm.habitat[mcm.habitat$Vegetation == "Desert",]
mcm.temperate.forest<-mcm.habitat[mcm.habitat$Vegetation == "Temperate forest",]
mcm.boreal.forest<-mcm.habitat[mcm.habitat$Vegetation== "Boreal forest",]


# use the unite function from package tidyr to merge the rank number from BS, Diet, and Locomotor
# codes here are retroprospect from the bootstrapping below

# resample.data.temp <- unite(input.data, "BS.Rank Diet.Rank Locomotor.Rank",c(BS.Rank, Diet.Rank, Locomotor.Rank), sep = "")

# resample function for speceis ecotypes most repeatedly
resample.func <- function (input.data, n.species)
{
  set.seed(1234)
  resample.data <- sample(input.data[,1], 1000, replace=T)
  table.data <- table(resample.data)
  order.data <- table.data[order(table.data, decreasing = T)][1:n.species]
  ## ordered data becaome data frame with two columns
  data.frame.resample<-data.frame(order.data)
  ## separate the BS, Diet, and Locomotor
  resample.data.type<-cbind(read.fwf(file=textConnection(as.character(data.frame.resample[,1])),
                                     widths = c(1,1,1), colClasses="character",
                                     col.names = c("BS.Rank","Diet.Rank","Locomotor.Rank")),
                            data.frame.resample[-1])
   return(resample.data.type)
}

# function for resampled ERich based on average No of species
resampled.ERich.func <- function (input.data)
  {
  resampled.ERich <- data.frame()
  resampled.data <- sample(input.data[,1], 6, replace=T)
  temp.resampled.data <- unique(resampled.data)
  resampled.ERich <- length(temp.resampled.data)
  resampled.ERich
}

# function to calculate mean of resampled number of speceis in each evironmental type
mean.species.func <- function (input.data)
{
  set.seed(1234)
  resampled.no.species <- sample(table(input.data[,1]), 1000, replace = T)
  mean.result <- mean(resampled.no.species)
  return(mean.result)
}

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
tc.type<-data.frame(mcm.tropical.climate$Eco.type)## tropical climate
ac.type<-data.frame(mcm.arid.climate$Eco.type) ## arid climate
ttc.type<-data.frame(mcm.temperate.climate$Eco.type) ## temperate climate
cc.type<-data.frame(mcm.cold.climate$Eco.type)## cold climate
trf.type<-data.frame(mcm.tropical.rainforest$Eco.type) ## tropical rianforest
tsf.type<-data.frame(mcm.tropical.seasonal.forest$Eco.type) ## tropical seasonal forest
sv.type<-data.frame(mcm.savanna$Eco.type)## savanna
gl.type<-data.frame(mcm.grassland$Eco.type)  ## grassland
sl.type<-data.frame(mcm.shrubland$Eco.type)## shrubland
dt.type<-data.frame(mcm.desert$Eco.type) ## desert
ttf.type<-data.frame(mcm.temperate.forest$Eco.type)## temperate forest
bf.type<-data.frame(mcm.boreal.forest$Eco.type)## boreal forest

# Two habiat
cl.resample.type<-resample.func(cl.type, n.species.close)
write.table(cl.resample.type[,-4],file = "sig.close_16Jan2018.csv", sep=",", row.names=F, col.names = F)
op.resample.type<-resample.func(op.type, n.species.open)
write.table(op.resample.type[,-4],file = "sig.open_16Jan2018.csv", sep=",", row.names=F, col.names = F)

# four climates
tc.resample.type<-resample.func(tc.type, n.species.tropical.climate)
write.table(tc.resample.type[,-4],file = "sig.tropical.climate_16Jan2018.csv", sep=",", row.names=F, col.names = F)
ac.resample.type<-resample.func(ac.type, n.species.arid.climate)
write.table(ac.resample.type[,-4],file = "sig.arid.climate_16Jan2018.csv", sep=",", row.names=F, col.names = F)
ttc.resample.type<-resample.func(ttc.type, n.species.temperate.climate)
write.table(ttc.resample.type[,-4],file = "sig.temperate.climate_16Jan2018.csv", sep=",", row.names=F, col.names = F)
cc.resample.type<-resample.func(cc.type, n.species.cold.climate)
write.table(cc.resample.type[,-4],file = "sig.cold.climate_16Jan2018.csv", sep=",", row.names=F, col.names = F)

# eight vegetations
trf.resample.type<-resample.func(trf.type, n.species.tropical.rainforest)
write.table(trf.resample.type[,-4],file = "sig.tropical.rainforest_16Jan2018.csv", sep=",", row.names=F, col.names = F)
tsf.resample.type<-resample.func(tsf.type, n.species.tropical.seasonal.forest)
write.table(tsf.resample.type[,-4],file = "sig.tropical.seasonal.forest_16Jan2018.csv", sep=",", row.names=F, col.names = F)
sv.resample.type<-resample.func(sv.type, n.species.savanna)
write.table(sv.resample.type[,-4],file = "sig.savanna_16Jan2018.csv", sep=",", row.names=F, col.names = F)
gl.resample.type<-resample.func(gl.type, n.species.grassland)
write.table(gl.resample.type[,-4],file = "sig.grassland_16Jan2018.csv", sep=",", row.names=F, col.names = F)
sl.resample.type<-resample.func(sl.type, n.species.shrubland)
write.table(sl.resample.type[,-4],file = "sig.shrubland_16Jan2018.csv", sep=",", row.names=F, col.names = F)
dt.resample.type<-resample.func(dt.type, n.species.desert)
write.table(dt.resample.type[,-4],file = "sig.desert_16Jan2018.csv", sep=",", row.names=F, col.names = F)
ttf.resample.type<-resample.func(ttf.type, n.species.temperate.forest)
write.table(ttf.resample.type[,-4],file = "sig.temperate.forest_16Jan2018.csv", sep=",", row.names=F, col.names = F)
bf.resample.type<-resample.func(bf.type, n.species.boreal.forest)
write.table(bf.resample.type[,-4],file = "sig.boreal.forest_16Jan2018.csv", sep=",", row.names=F, col.names = F)

##
sig.mastersheet <- read.csv(file = "Sig_MaterSheet_Resampled_16Jan2018.csv", sep = ",")
sig.mastersheet <- sig.mastersheet[,c(1,5)]

write.table(table(sig.mastersheet), file="Sig_MaterSheet_Summary_Resampled_16Jan2018.csv", sep=",")

########################################################################
### Bootstrapping for designated numbers of species with a community ###
########################################################################

###############################################################################
### Don't run it again because the resampling data has been saved
### If you accidentally run the resamppling code, the plots and distribution of
### may be shifted a little bit. However, all plots are needed to be revised.
###############################################################################

##############  Bootstrapping function

BS.Species.func  <- function (input.data, n.species)
{
  mean.eds<-data.frame()
  set.seed(1234)
  for (t in 1:1000)
  {
    eds.community<-data.frame()
    data.temp<-input.data[sample(nrow(input.data), n.species, replace=T),]
    for (u in 1:(n.species-1))
    {
      for (v in u:(n.species-1))
      {
        results<-abs(data.temp[u,]-data.temp[v+1,])
        eds.temp<-rowSums(results)
        eds.community<-rbind(eds.community,c(eds.temp))
      }
    }
    mean.eds<-rbind(mean.eds,c(mean(eds.community[,1])))
  }
  return(mean.eds)
}

#### alternaitive function
# BS.Species.func  <- function (input.data, n.species)
# {
#   eds.community<-data.frame()
#   for (t in 1:1000)
#   {
#     data.temp<-input.data[sample(nrow(input.data), n.species, replace=F),]
#     for (u in 1:(n.species-1))
#     {
#       for (v in u:(n.species-1))
#       {
#         results<-abs(data.temp[u,]-data.temp[v+1,])
#         eds.temp<-rowSums(results)
#         eds.community<-rbind(eds.community,c(eds.temp))
#       }
#     }
#   }
#   return(eds.community)
# }

###########################################
## 5 Species in each simulated communities
###########################################

##################################
### 2 habitats
## close habitat
start<-Sys.time()
cl.results<-BS.Species.func(mcm.close[,3:5], 5)
print(Sys.time()-start)

## open
start<-Sys.time()
op.results<-BS.Species.func(mcm.open[,3:5], 5)
print(Sys.time()-start)


##################################
### 4 different climates

## Tropical climate
start<-Sys.time()
tc.results<-BS.Species.func(mcm.tropical.climate[,3:5], 5)
print(Sys.time()-start)

## Arid climate
start<-Sys.time()
ac.results<-BS.Species.func(mcm.arid.climate[,3:5], 5)
print(Sys.time()-start)

## Temperate climate
start<-Sys.time()
ttc.results<-BS.Species.func(mcm.temperate.climate[,3:5], 5)
print(Sys.time()-start)

## Cold climate
start<-Sys.time()
cc.results<-BS.Species.func(mcm.cold.climate[,3:5], 5)
print(Sys.time()-start)

### Eight Vegetation Structure
### Tropical rainforest
start<-Sys.time()
trf.results<-BS.Species.func(mcm.tropical.rainforest[,3:5], 5)
print(Sys.time()-start)

### Tropical seasonal forest
start<-Sys.time()
tsf.results<-BS.Species.func(mcm.tropical.seasonal.forest[,3:5], 5)
print(Sys.time()-start)

## Savanna
start<-Sys.time()
sv.results<-BS.Species.func(mcm.savanna[,3:5], 5)
print(Sys.time()-start)

## grassland
start<-Sys.time()
gl.results<-BS.Species.func(mcm.grassland[,3:5], 5)
print(Sys.time()-start)

## Shrubland
start<-Sys.time()
sl.results<-BS.Species.func(mcm.shrubland[,3:5], 5)
print(Sys.time()-start)

## Desert
start<-Sys.time()
dt.results<-BS.Species.func(mcm.desert[,3:5], 5)
print(Sys.time()-start)

## Temperate forest
start<-Sys.time()
ttf.results<-BS.Species.func(mcm.temperate.forest[,3:5], 5)
print(Sys.time()-start)

## Boreal forest
start<-Sys.time()
bf.results<-BS.Species.func(mcm.boreal.forest[,3:5], 5)
print(Sys.time()-start)

## Mean and SD of all bootstrapped values
bootstrapped.results.all.5p<-matrix(c(cl.results[,1], op.results[,1],
                                   tc.results[,1], ac.results[,1], ttc.results[,1], cc.results[,1],
                                   trf.results[,1], tsf.results[,1], sv.results[,1], gl.results[,1],
                                   sl.results[,1], dt.results[,1], ttf.results[,1], bf.results[,1]),
                                 nrow = 1000, ncol = 14, byrow = F)

colnames(bootstrapped.results.all.5p)<-c("Close","Open",
                                      "Tropical","Arid", "Temperate", "Cold",
                                      "Tropical rainforest", "Tropical seasonal forest", "Savanna","Grassland",
                                      "Shrubland", "Desert", "Temperate forest", "Boreal forest")

mean.bootstrapped.5p<-apply(bootstrapped.results.all.5p, 2, mean)
sd.bootstrapped.5p<-apply(bootstrapped.results.all.5p, 2, sd)
qt.bootstrapped.5p<-apply(bootstrapped.results.all.5p, 2, quantile, prob=c(0.025, 0.975))

## manage the data
bootstropped.results.stat.5p <- t(rbind(t(mean.bootstrapped.5p), t(sd.bootstrapped.5p)))
colnames(bootstropped.results.stat.5p) <- c("Mean", "sd")
bootstropped.results.stat.5p <- cbind(bootstropped.results.stat.5p, t(qt.bootstrapped.5p))

write.csv(bootstropped.results.stat.5p, file = "Bootstrapped_statistic_5species_16Jan2018.csv")

##########################
### Pair-wise t test  ####
##########################
habitat.pair <- c("open vs close")
climate.pair <- c ("tropical vs arid", "tropical vs temperate", "tropical vs cold",
                   "arid vs temperate", "arid vs cold",
                   "temperate vs cold")
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

# ## p value funrctions
# ## p-value
# func.p.value <- function (bootstrap.a, bootstrap.b, obs.diff)
# {
#   diff <- rep(NA, 1000)
#   for (i in 1000) {
#     diff[i] <- bootstrap.a[i]-bootstrap.b[i]
#   }
#     ## two-tail p-value; -0.15 = mean(dep.close) - mean(dep.open)
#   p.value <- cat("p<", length(diff[diff>=obs.diff|diff <= obs.diff])/length(diff), "\n")
#   return(p.value)
# }


##### two habitats
habitat.results.test.5p <- data.frame()
habitat.results.test.5p<-t.test(bootstrapped.results.all.5p[,1], bootstrapped.results.all.5p[,2])
habitat.results.test.5p<-rbind(c(habitat.results.test.5p$statistic,
                              habitat.results.test.5p$parameter,
                              habitat.results.test.5p$p.value))
colnames(habitat.results.test.5p)<-c("t","df","P value")
rownames(habitat.results.test.5p) <- habitat.pair

##### four climate
climate.results.test.5p<-data.frame()

for (z in 3:5)
{
  for (h in z:5) {
    temp.result.test.5p<-t.test(bootstrapped.results.all.5p[,z], bootstrapped.results.all.5p[, h+1])
    climate.results.test.5p<-rbind(climate.results.test.5p, c(temp.result.test.5p$statistic, temp.result.test.5p$parameter, temp.result.test.5p$p.value))
  }
}

colnames(climate.results.test.5p)<-c("t","df","P value")
rownames(climate.results.test.5p) <- climate.pair

#### seven vegetations
vegetation.results.test.5p<-data.frame()

for (w in 7:13)
{
  for (q in w:13) {
    temp.result.test.5p<-t.test(bootstrapped.results.all.5p[,w], bootstrapped.results.all.5p[, q+1])
    vegetation.results.test.5p<-rbind(vegetation.results.test.5p,c(temp.result.test.5p$statistic,temp.result.test.5p$parameter,temp.result.test.5p$p.value))
  }
}

colnames(vegetation.results.test.5p)<-c("t","df","P value")
rownames(vegetation.results.test.5p) <- vegetation.pair

## combine all results
t.test.results.all.5species <- rbind(habitat.results.test.5p,
                                     climate.results.test.5p,
                                     vegetation.results.test.5p)

write.csv(t.test.results.all.5species, file = "Bootstrapped_statistic_table_5species_16Jan2018.csv")

# ##### t.test with JLS (5-species community)
# t.test.results.JLS<-data.frame()
#
# for (i in 1:14) {
#   t.test.results.JLS.temp<-t.test(bootstrapped.results.all[,i], alternative="two.sided", mu = mean(eds.jls[,1]))
#   t.test.results.JLS<-rbind(t.test.results.JLS, c(t.test.results.JLS.temp$statistic,
#                                                   t.test.results.JLS.temp$parameter,
#                                                   t.test.results.JLS.temp$p.value))
# }
#
# colnames(t.test.results.JLS)<-c("t","df","P value")
# # write.csv(t.test.results.JLS, file = "Bootstrapped statistic table JLS 5_species.csv")

## histgram plots
par(mfcol=c(8,3),
    oma=c(1,1,0,1),
    mar=c(2,2,1,0),
    mgp=c(1,0.6,0),
    tck=-0.02)

# 2 habitats
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space

h<-hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey30")
xfit<-seq(min(cl.results[,1]), max(cl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(cl.results[,1]),sd=sd(cl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(cl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(cl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="black")
xfit<-seq(min(op.results[,1]), max(op.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(op.results[,1]),sd=sd(op.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(op.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(op.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space

# 4 climates
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space

h<-hist(tc.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkorange")
xfit<-seq(min(tc.results[,1]), max(tc.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(tc.results[,1]),sd=sd(tc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(tc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(tc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)


h<-hist(ac.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey20")
xfit<-seq(min(ac.results[,1]), max(ac.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(ac.results[,1]),sd=sd(ac.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ac.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ac.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(ttc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
xfit<-seq(min(ttc.results[,1]), max(ttc.results[,1]), length=40)
yfit<-dnorm(xfit, mean=mean(ttc.results[,1]),sd=sd(ttc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ttc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ttc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),
        xlab="Bootstrapped ecological disparity", ylab="", main="", col="dodgerblue3")
xfit<-seq(min(cc.results[,1]), max(cc.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(cc.results[,1]),sd=sd(cc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(cc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(cc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space

## 8 vegetations
h<-hist(trf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="red")
xfit<-seq(min(trf.results[,1]), max(trf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(trf.results[,1]),sd=sd(trf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(trf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(trf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(tsf.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="orangered")
xfit<-seq(min(tsf.results[,1]), max(tsf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(tsf.results[,1]),sd=sd(tsf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(tsf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(tsf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(sv.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="firebrick1")
xfit<-seq(min(sv.results[,1]), max(sv.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(sv.results[,1]),sd=sd(sv.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(sv.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(sv.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(gl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="goldenrod")
xfit<-seq(min(gl.results[,1]), max(gl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(gl.results[,1]),sd=sd(gl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(gl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(gl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(sl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkgoldenrod1")
xfit<-seq(min(sl.results[,1]), max(sl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(sl.results[,1]),sd=sd(sl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(sl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(sl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(dt.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col= "khaki")
xfit<-seq(min(dt.results[,1]), max(dt.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(dt.results[,1]),sd=sd(dt.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(dt.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(dt.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(ttf.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
xfit<-seq(min(ttf.results[,1]), max(ttf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(ttf.results[,1]),sd=sd(ttf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ttf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ttf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(bf.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkviolet")
xfit<-seq(min(bf.results[,1]), max(bf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(bf.results[,1]),sd=sd(bf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(bf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(bf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)


###########################################
## 6 Species in each simulated communities
################################

##################################
### 2 habitats
## close habitat
start<-Sys.time()
cl.results<-BS.Species.func(mcm.close[,3:5], 6)
print(Sys.time()-start)

## open
start<-Sys.time()
op.results<-BS.Species.func(mcm.open[,3:5], 6)
print(Sys.time()-start)


##################################
### 4 different climates

## Tropical climate
start<-Sys.time()
tc.results<-BS.Species.func(mcm.tropical.climate[,3:5], 6)
print(Sys.time()-start)

## Arid climate
start<-Sys.time()
ac.results<-BS.Species.func(mcm.arid.climate[,3:5], 6)
print(Sys.time()-start)

## Temperate climate
start<-Sys.time()
ttc.results<-BS.Species.func(mcm.temperate.climate[,3:5], 6)
print(Sys.time()-start)

## Cold climate
start<-Sys.time()
cc.results<-BS.Species.func(mcm.cold.climate[,3:5], 6)
print(Sys.time()-start)

### Seven Vegetation Structure
### Tropical rainforest
start<-Sys.time()
trf.results<-BS.Species.func(mcm.tropical.rainforest[,3:5], 6)
print(Sys.time()-start)

### Tropical seasonal forest
start<-Sys.time()
tsf.results<-BS.Species.func(mcm.tropical.seasonal.forest[,3:5], 6)
print(Sys.time()-start)

## savanna
start<-Sys.time()
sv.results<-BS.Species.func(mcm.savanna[,3:5], 6)
print(Sys.time()-start)

## grassland
start<-Sys.time()
gl.results<-BS.Species.func(mcm.grassland[,3:5], 6)
print(Sys.time()-start)

## Shrubland
start<-Sys.time()
sl.results<-BS.Species.func(mcm.shrubland[,3:5], 6)
print(Sys.time()-start)

## Desert
start<-Sys.time()
dt.results<-BS.Species.func(mcm.desert[,3:5], 6)
print(Sys.time()-start)

## Temperate forest
start<-Sys.time()
ttf.results<-BS.Species.func(mcm.temperate.forest[,3:5], 6)
print(Sys.time()-start)

## Boreal forest
start<-Sys.time()
bf.results<-BS.Species.func(mcm.boreal.forest[,3:5], 6)
print(Sys.time()-start)

## Mean and SD of all bootstrapped values
bootstrapped.results.all.6p<-matrix(c(cl.results[,1], op.results[,1],
                                   tc.results[,1], ac.results[,1], ttc.results[,1], cc.results[,1],
                                   trf.results[,1], tsf.results[,1], sv.results[,1], gl.results[,1],
                                   sl.results[,1], dt.results[,1], ttf.results[,1], bf.results[,1]),
                                 nrow = 1000, ncol = 14, byrow = F)

colnames(bootstrapped.results.all.6p)<-c("Close","Open",
                                      "Tropical","Arid", "Temperate", "Cold",
                                      "Tropical rainforest", "Tropical seasonal forest", "Savanna","Grassland",
                                      "Shrubland", "Desert", "Temperate forest", "Boreal forest")

mean.bootstrapped.6p<-apply(bootstrapped.results.all.6p, 2, mean)
sd.bootstrapped.6p<-apply(bootstrapped.results.all.6p, 2, sd)
qt.bootstrapped.6p<-apply(bootstrapped.results.all.6p, 2, quantile, prob=c(0.025, 0.975))

## manage the data
bootstropped.results.stat.6p <- t(rbind(t(mean.bootstrapped.6p), t(sd.bootstrapped.6p)))
colnames(bootstropped.results.stat.6p) <- c("Mean", "sd")
bootstropped.results.stat.6p <- cbind(bootstropped.results.stat.6p, t(qt.bootstrapped.6p))

write.csv(bootstropped.results.stat.6p, file = "Bootstrapped_statistic_6species_16Jan2018.csv")

##########################
### Pair-wise t test  ####
##########################
##### two habitats
habitat.results.test.6p <- data.frame()
habitat.results.test.6p<-t.test(bootstrapped.results.all.6p[,1], bootstrapped.results.all.6p[,2])
habitat.results.test.6p<-rbind(c(habitat.results.test.6p$statistic,
                              habitat.results.test.6p$parameter,
                              habitat.results.test.6p$p.value))
colnames(habitat.results.test.6p)<-c("t","df","P value")
rownames(habitat.results.test.6p) <- habitat.pair

##### four climate
climate.results.test.6p<-data.frame()

for (z in 3:5)
{
  for (h in z:5) {
    temp.result.test.6p<-t.test(bootstrapped.results.all.6p[,z], bootstrapped.results.all.6p[, h+1])
    climate.results.test.6p<-rbind(climate.results.test.6p,c(temp.result.test.6p$statistic,temp.result.test.6p$parameter,temp.result.test.6p$p.value))
  }
}

colnames(climate.results.test.6p)<-c("t","df","P value")
rownames(climate.results.test.6p) <- climate.pair

## p value

#### eight vegetations
vegetation.results.test.6p<-data.frame()

for (w in 7:13)
{
  for (q in w:13) {
    temp.result.test.6p<-t.test(bootstrapped.results.all.6p[,w], bootstrapped.results.all.6p[, q+1])
    vegetation.results.test.6p<-rbind(vegetation.results.test.6p,c(temp.result.test.6p$statistic,temp.result.test.6p$parameter,temp.result.test.6p$p.value))
  }
}

colnames(vegetation.results.test.6p)<-c("t","df","P value")
rownames(vegetation.results.test.6p) <- vegetation.pair

## combine all results
t.test.results.all.6species <- rbind(habitat.results.test.6p,
                                     climate.results.test.6p,
                                     vegetation.results.test.6p)

write.csv(t.test.results.all.6species, file = "Bootstrapped_statistic_table_6species_16Jan2018.csv")

# ##### t.test with DWZZ (6-species community)
# t.test.results.DWZZ<-data.frame()
#
# for (i in 1:14) {
#   t.test.results.DWZZ.temp<-t.test(bootstrapped.results.all[,i], alternative="two.sided", mu = mean(eds.dwzz[,1]))
#   t.test.results.DWZZ<-rbind(t.test.results.DWZZ, c(t.test.results.DWZZ.temp$statistic,
#                                                   t.test.results.DWZZ.temp$parameter,
#                                                   t.test.results.DWZZ.temp$p.value))
# }
#
# colnames(t.test.results.DWZZ)<-c("t","df","P value")
# # write.csv(t.test.results.DWZZ, file = "Bootstrapped statistic table DWZZ 5_species.csv")

## histgram plots
par(mfcol=c(8,3),
    oma=c(1,1,0,1),
    mar=c(2,2,1,0),
    mgp=c(1,0.6,0),
    tck=-0.02)

# 2 habitats
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space

h<-hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey30")
xfit<-seq(min(cl.results[,1]), max(cl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(cl.results[,1]),sd=sd(cl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(cl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(cl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(op.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="black")
xfit<-seq(min(op.results[,1]), max(op.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(op.results[,1]),sd=sd(op.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(op.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(op.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space

# 4 climates
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space

h<-hist(tc.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkorange")
xfit<-seq(min(tc.results[,1]), max(tc.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(tc.results[,1]),sd=sd(tc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(tc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(tc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)


h<-hist(ac.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey20")
xfit<-seq(min(ac.results[,1]), max(ac.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(ac.results[,1]),sd=sd(ac.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ac.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ac.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(ttc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
xfit<-seq(min(ttc.results[,1]), max(ttc.results[,1]), length=40)
yfit<-dnorm(xfit, mean=mean(ttc.results[,1]),sd=sd(ttc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ttc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ttc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(cc.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120),
        xlab="Bootstrapped ecological disparity", ylab="", main="", col="dodgerblue3")
xfit<-seq(min(cc.results[,1]), max(cc.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(cc.results[,1]),sd=sd(cc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(cc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(cc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space

## 8 vegetations
h<-hist(trf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="red")
xfit<-seq(min(trf.results[,1]), max(trf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(trf.results[,1]),sd=sd(trf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(trf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(trf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(tsf.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="orangered")
xfit<-seq(min(tsf.results[,1]), max(tsf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(tsf.results[,1]),sd=sd(tsf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(tsf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(tsf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(sv.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F, col="firebrick1")
xfit<-seq(min(sv.results[,1]), max(sv.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(sv.results[,1]),sd=sd(sv.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(sv.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(sv.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(gl.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="goldenrod")
xfit<-seq(min(gl.results[,1]), max(gl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(gl.results[,1]),sd=sd(gl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(gl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(gl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(sl.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkgoldenrod1")
xfit<-seq(min(sl.results[,1]), max(sl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(sl.results[,1]),sd=sd(sl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(sl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(sl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(dt.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col= "khaki")
xfit<-seq(min(dt.results[,1]), max(dt.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(dt.results[,1]),sd=sd(dt.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(dt.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(dt.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(ttf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
xfit<-seq(min(ttf.results[,1]), max(ttf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(ttf.results[,1]),sd=sd(ttf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ttf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ttf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(bf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkviolet")
xfit<-seq(min(bf.results[,1]), max(bf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(bf.results[,1]),sd=sd(bf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(bf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(bf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)


###########################################
## 10 Species in each simulated communities
################################
### 2 habitats
## close habitat
start<-Sys.time()
cl.results<-BS.Species.func(mcm.close[,3:5], 10)
print(Sys.time()-start)

## open
start<-Sys.time()
op.results<-BS.Species.func(mcm.open[,3:5], 10)
print(Sys.time()-start)


##################################
### 4 different climates

## Tropical climate
start<-Sys.time()
tc.results<-BS.Species.func(mcm.tropical.climate[,3:5], 10)
print(Sys.time()-start)

## Arid climate
start<-Sys.time()
ac.results<-BS.Species.func(mcm.arid.climate[,3:5], 10)
print(Sys.time()-start)

## Temperate climate
start<-Sys.time()
ttc.results<-BS.Species.func(mcm.temperate.climate[,3:5], 10)
print(Sys.time()-start)

## Cold climate
start<-Sys.time()
cc.results<-BS.Species.func(mcm.cold.climate[,3:5], 10)
print(Sys.time()-start)

### Seven Vegetation Structure
### Tropical rainforest
start<-Sys.time()
trf.results<-BS.Species.func(mcm.tropical.rainforest[,3:5], 10)
print(Sys.time()-start)

### Tropical seasonal forest
start<-Sys.time()
tsf.results<-BS.Species.func(mcm.tropical.seasonal.forest[,3:5], 10)
print(Sys.time()-start)

## savanna
start<-Sys.time()
sv.results<-BS.Species.func(mcm.savanna[,3:5], 10)
print(Sys.time()-start)

## grassland
start<-Sys.time()
gl.results<-BS.Species.func(mcm.grassland[,3:5], 10)
print(Sys.time()-start)

## Shrubland
start<-Sys.time()
sl.results<-BS.Species.func(mcm.shrubland[,3:5], 10)
print(Sys.time()-start)

## Desert
start<-Sys.time()
dt.results<-BS.Species.func(mcm.desert[,3:5], 10)
print(Sys.time()-start)

## Temperate forest
start<-Sys.time()
ttf.results<-BS.Species.func(mcm.temperate.forest[,3:5], 10)
print(Sys.time()-start)

## Boreal forest
start<-Sys.time()
bf.results<-BS.Species.func(mcm.boreal.forest[,3:5], 10)
print(Sys.time()-start)

## Mean and SD of all bootstrapped values
bootstrapped.results.all.10p<-matrix(c(cl.results[,1], op.results[,1],
                                   tc.results[,1], ac.results[,1], ttc.results[,1], cc.results[,1],
                                   trf.results[,1], tsf.results[,1], sv.results[,1], gl.results[,1],
                                   sl.results[,1], dt.results[,1], ttf.results[,1], bf.results[,1]),
                                 nrow = 1000, ncol = 14, byrow = F)

colnames(bootstrapped.results.all.10p)<-c("Close","Open",
                                      "Tropical","Arid", "Temperate", "Cold",
                                      "Tropical rainforest", "Tropical seasonal forest", "Savanna","Grassland",
                                      "Shrubland", "Desert", "Temperate forest", "Boreal forest")

mean.bootstrapped.10p<-apply(bootstrapped.results.all.10p, 2, mean)
sd.bootstrapped.10p<-apply(bootstrapped.results.all.10p, 2, sd)
qt.bootstrapped.10p<-apply(bootstrapped.results.all.10p, 2, quantile, prob=c(0.025, 0.975))

## manage the data
bootstropped.results.stat.10p <- t(rbind(t(mean.bootstrapped.10p), t(sd.bootstrapped.10p)))
colnames(bootstropped.results.stat.10p) <- c("Mean", "sd")
bootstropped.results.stat.10p <- cbind(bootstropped.results.stat.10p, t(qt.bootstrapped.10p))

write.csv(bootstropped.results.stat.10p, file = "Bootstrapped_statistic_10species_16Jan2018.csv")

##########################
### Pair-wise t test  ####
##########################
##### two habitats
habitat.results.test.10p <- data.frame()
habitat.results.test.10p <-t.test(bootstrapped.results.all.10p[,1], bootstrapped.results.all.10p[,2])
habitat.results.test.10p <-rbind(c(habitat.results.test.10p$statistic,
                              habitat.results.test.10p$parameter,
                              habitat.results.test.10p$p.value))
colnames(habitat.results.test.10p)<-c("t","df","P value")
rownames(habitat.results.test.10p) <- habitat.pair

##### four climate
climate.results.test.10p<-data.frame()

for (z in 3:5)
{
  for (h in z:5) {
    temp.result.test.10p<-t.test(bootstrapped.results.all.10p[,z], bootstrapped.results.all.10p[, h+1])
    climate.results.test.10p<-rbind(climate.results.test.10p,c(temp.result.test.10p$statistic,temp.result.test.10p$parameter,temp.result.test.10p$p.value))
  }
}

colnames(climate.results.test.10p)<-c("t","df","P value")
rownames(climate.results.test.10p) <- climate.pair

#### eight vegetations
vegetation.results.test.10p<-data.frame()

for (w in 7:13)
{
  for (q in w:13) {
    temp.result.test.10p<-t.test(bootstrapped.results.all.10p[,w], bootstrapped.results.all.10p[, q+1])
    vegetation.results.test.10p<-rbind(vegetation.results.test.10p,c(temp.result.test.10p$statistic,temp.result.test.10p$parameter,temp.result.test.10p$p.value))
  }
}

colnames(vegetation.results.test.10p)<-c("t","df","P value")
rownames(vegetation.results.test.10p) <- vegetation.pair

## combine all results
t.test.results.all.10species <- rbind(habitat.results.test.10p,
                                     climate.results.test.10p,
                                     vegetation.results.test.10p)

write.csv(t.test.results.all.10species, file = "Bootstrapped_statistic_table_10species_16Jan2018.csv")

# ##### t.test with LJT_JSG (11-species community)
# t.test.results.LJT_JSG<-data.frame()
#
# for (i in 1:14) {
#   t.test.results.LJT_JSG.temp<-t.test(bootstrapped.results.all[,i], alternative="two.sided", mu = mean(eds.ljt.jsg[,1]))
#   t.test.results.LJT_JSG<-rbind(t.test.results.LJT_JSG, c(t.test.results.LJT_JSG.temp$statistic,
#                                                   t.test.results.LJT_JSG.temp$parameter,
#                                                   t.test.results.LJT_JSG.temp$p.value))
# }
#
# colnames(t.test.results.LJT_JSG)<-c("t","df","P value")
# write.csv(t.test.results.LJT_JSG, file = "Bootstrapped statistic table LJT_JSG 7_species.csv")

## histgram plots
par(mfcol=c(8,3),
    oma=c(1,1,0,1),
    mar=c(2,2,1,0),
    mgp=c(1,0.6,0),
    tck=-0.02)

# 2 habitats
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space

h<-hist(cl.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey30")
xfit<-seq(min(cl.results[,1]), max(cl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(cl.results[,1]),sd=sd(cl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(cl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(cl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(op.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,120), ann=F, col="black")
xfit<-seq(min(op.results[,1]), max(op.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(op.results[,1]),sd=sd(op.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(op.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(op.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space

# 4 climates
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space

h<-hist(tc.results[,1], breaks = 25, xlim = c(0,8), ylim = c(0,150), ann=F, col="darkorange")
xfit<-seq(min(tc.results[,1]), max(tc.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(tc.results[,1]),sd=sd(tc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(tc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(tc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)


h<-hist(ac.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey20")
xfit<-seq(min(ac.results[,1]), max(ac.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(ac.results[,1]),sd=sd(ac.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ac.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ac.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(ttc.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
xfit<-seq(min(ttc.results[,1]), max(ttc.results[,1]), length=40)
yfit<-dnorm(xfit, mean=mean(ttc.results[,1]),sd=sd(ttc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ttc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ttc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(cc.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,150),
        xlab="Bootstrapped ecological disparity", ylab="", main="", col="dodgerblue3")
xfit<-seq(min(cc.results[,1]), max(cc.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(cc.results[,1]),sd=sd(cc.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(cc.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(cc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space

## 8 vegetations
h<-hist(trf.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,150), ann=F, col="red")
xfit<-seq(min(trf.results[,1]), max(trf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(trf.results[,1]),sd=sd(trf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(trf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(trf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(tsf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,150), ann=F, col="orangered")
xfit<-seq(min(tsf.results[,1]), max(tsf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(tsf.results[,1]),sd=sd(tsf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(tsf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(tsf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(sv.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="firebrick1")
xfit<-seq(min(sv.results[,1]), max(sv.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(sv.results[,1]),sd=sd(sv.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(sv.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(sv.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(gl.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,180), ann=F, col="goldenrod")
xfit<-seq(min(gl.results[,1]), max(gl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(gl.results[,1]),sd=sd(gl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(gl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(gl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(sl.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,150), ann=F, col="darkgoldenrod1")
xfit<-seq(min(sl.results[,1]), max(sl.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(sl.results[,1]),sd=sd(sl.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(sl.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(sl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(dt.results[,1], breaks = 15, xlim = c(0,8), ylim = c(0,200), ann=F, col= "khaki")
xfit<-seq(min(dt.results[,1]), max(dt.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(dt.results[,1]),sd=sd(dt.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(dt.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(dt.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(ttf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,150), ann=F, col="green3")
xfit<-seq(min(ttf.results[,1]), max(ttf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(ttf.results[,1]),sd=sd(ttf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(ttf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(ttf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

h<-hist(bf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,150), ann=F, col="darkviolet")
xfit<-seq(min(bf.results[,1]), max(bf.results[,1]), length=40)
yfit<-dnorm(xfit,mean=mean(bf.results[,1]),sd=sd(bf.results[,1]))
yfit <- yfit*diff(h$mids[1:2])*length(bf.results[,1])
lines(xfit, yfit, col="blue", lwd=1)
abline(v=quantile(bf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)

# ############################################################################
# ##### t.test with JLS (5-species community)
# t.test.results.JLS<-data.frame()
#
# ##### t.test with DWZZ (6-species community)
# t.test.results.DWZZ<-data.frame()
#
# ##### t.test with TJS (7-species community)
# t.test.results.TJS<-data.frame()
#
# ##### t.test with LJT-JSG (11-species community)
# t.test.results.LJT_JSG<-data.frame()

ChMeMaCo<-read.csv(file="ChinaMesozoicMammalCommunity_28June2018.csv", header = T, sep=",")

### Resample the same number of species of Mesozoic mammalian communtiies with regard to the Messel community

Messel.community <- ChMeMaCo[ChMeMaCo$Community=="MESL",]
Messel.community$Eco.type <- with(Messel.community, paste0(BSRank, DietRank, LocomotorRank))
Messel.ecotype <- Messel.community$Eco.type

resample.func <- function (input.data, n.species, r)
{
  resample.data <- mat.or.vec(nr = r, nc = n.species )
  for (i in 1:r)
    {
    resample.data[i,] <- sample(input.data, n.species, replace=T)
  }
  return(resample.data)
}


### DWZZ
Messel.DWZZ.Resample.ecotype <- resample.func(Messel.ecotype, 6, 10000)
Messel.DWZZ.Resample.ecotype.sorted <- sort(table(Messel.DWZZ.Resample.ecotype), decreasing = T)
Messel.DWZZ.Resample.ecotype.sorted.6species <- Messel.DWZZ.Resample.ecotype.sorted[1:6]
write.csv(Messel.DWZZ.Resample.ecotype.sorted.6species, file = "Resampled_DWZZ_by_Messel_29June2018.csv")

### LJT-JSG and TJS would be the same becaues they have the same number of speccies in each community
Messel.LJT_JLS.Resample.ecotype <- resample.func(Messel.ecotype, 10, 10000)
Messel.LJT_JLS.Resample.ecotype.sorted <- sort(table(Messel.LJT_JLS.Resample.ecotype), decreasing = T)
Messel.LJT_JLS.Resample.ecotype.sorted.10species <- Messel.LJT_JLS.Resample.ecotype.sorted[1:10]
write.csv(Messel.LJT_JLS.Resample.ecotype.sorted.10species, file = "Resampled_LJT-JSG_TJS_by_Messel_29June2018.csv")

### JLS
Messel.JSL.Resample.ecotype <- resample.func(Messel.ecotype, 5, 10000)
Messel.JSL.Resample.ecotype.sorted <- sort(table(Messel.JSL.Resample.ecotype), decreasing = T)
Messel.JSL.Resample.ecotype.sorted.5species <- Messel.JSL.Resample.ecotype.sorted[1:5]
write.csv(Messel.JSL.Resample.ecotype.sorted.5species, file = "Resampled_JSL_by_Messel_29June2018.csv")

### resampled the extinct mammalian communiites to look at the EDisp and ERcih distribution

extinct.community <- unique(ChMeMaCo$Community)
mean.resampled.edisp.extinct <- data.frame(matrix(NA, nrow = 5, ncol = 1000))
edisp.extinct <- data.frame()

for (dj in 1:length(extinct.community)) {
  temp.extinct.community.data <- data.frame()
  temp.extinct.community.data <- ChMeMaCo[ChMeMaCo$Community==extinct.community[dj],]
  n.extinct.species <- dim(temp.extinct.community.data)[1]
  species.position <- temp.extinct.community.data$SpeciesNo

  for (aj in 1:1000) {

    resampled.extinct.community.data <- data.frame()

    resampled.species.position<- sample(species.position, n.extinct.species, replace = T)

    for (zj in 1:n.extinct.species) {
      temp.resampled.extinct.community.data <- temp.extinct.community.data[temp.extinct.community.data$SpeciesNo==resampled.species.position[zj], ]
      resampled.extinct.community.data <- rbind(resampled.extinct.community.data, c(temp.resampled.extinct.community.data))
    }
    for (i in 1:(n.extinct.species-1)) {
      for (j in i:(n.extinct.species-1)) {
        results.extinct.abs<-abs(resampled.extinct.community.data[,4:6][i,]-resampled.extinct.community.data[,4:6][j+1,])
        edisp.extinct.temp<-rowSums(results.extinct.abs)
        edisp.extinct<-rbind(edisp.extinct,c(edsip.extinct.temp))
      }
    }
    mean.resampled.edisp.extinct[dj, aj]<-c(mean(edisp.extinct[,1]))
  }
}


## Backup codes for 9 species

# ###########################################
# ## 9 Species in each simulated communities
# ################################
#
# ### 2 habitats
# ## close habitat
# start<-Sys.time()
# cl.results<-BS.Species.func(mcm.close[,3:5], 9)
# print(Sys.time()-start)
#
# ## open
# start<-Sys.time()
# op.results<-BS.Species.func(mcm.open[,3:5], 9)
# print(Sys.time()-start)
#
#
# ##################################
# ### 4 different climates
#
# ## Tropical climate
# start<-Sys.time()
# tc.results<-BS.Species.func(mcm.tropical.climate[,3:5], 9)
# print(Sys.time()-start)
#
# ## Arid climate
# start<-Sys.time()
# ac.results<-BS.Species.func(mcm.arid.climate[,3:5], 9)
# print(Sys.time()-start)
#
# ## Temperate climate
# start<-Sys.time()
# ttc.results<-BS.Species.func(mcm.temperate.climate[,3:5], 9)
# print(Sys.time()-start)
#
# ## Cold climate
# start<-Sys.time()
# cc.results<-BS.Species.func(mcm.cold.climate[,3:5], 9)
# print(Sys.time()-start)
#
# ### Seven Vegetation Structure
# ### Tropical rainforest
# start<-Sys.time()
# trf.results<-BS.Species.func(mcm.tropical.rainforest[,3:5], 9)
# print(Sys.time()-start)
#
# ### Tropical seasonal forest
# start<-Sys.time()
# tsf.results<-BS.Species.func(mcm.tropical.seasonal.forest[,3:5], 9)
# print(Sys.time()-start)
#
# ## savanna
# start<-Sys.time()
# sv.results<-BS.Species.func(mcm.savanna[,3:5], 9)
# print(Sys.time()-start)
#
# ## grassland
# start<-Sys.time()
# gl.results<-BS.Species.func(mcm.grassland[,3:5], 9)
# print(Sys.time()-start)
#
# ## Shrubland
# start<-Sys.time()
# sl.results<-BS.Species.func(mcm.shrubland[,3:5], 9)
# print(Sys.time()-start)
#
# ## Desert
# start<-Sys.time()
# dt.results<-BS.Species.func(mcm.desert[,3:5], 9)
# print(Sys.time()-start)
#
# ## Temperate forest
# start<-Sys.time()
# ttf.results<-BS.Species.func(mcm.temperate.forest[,3:5], 9)
# print(Sys.time()-start)
#
# ## Boreal forest
# start<-Sys.time()
# bf.results<-BS.Species.func(mcm.boreal.forest[,3:5], 9)
# print(Sys.time()-start)
#
# ## Mean and SD of all bootstrapped values
# bootstrapped.results.all.9p<-matrix(c(cl.results[,1], op.results[,1],
#                                    tc.results[,1], ac.results[,1], ttc.results[,1], cc.results[,1],
#                                    trf.results[,1], tsf.results[,1], sv.results[,1], gl.results[,1], sl.results[,1], dt.results[,1], ttf.results[,1], bf.results[,1]),
#                                  nrow = 1000, ncol = 14, byrow = F)
#
# colnames(bootstrapped.results.all.9p)<-c("Close","Open",
#                                       "Tropical","Arid", "Temperate", "Cold",
#                                       "Tropical rainforest", "Tropical seasonal forest", "Savanna","Grassland", "Shrubland", "Desert", "Temperate forest", "Boreal forest")
#
# mean.bootstrapped.9p<-apply(bootstrapped.results.all.9p, 2, mean)
# sd.bootstrapped.9p<-apply(bootstrapped.results.all.9p, 2, sd)
# qt.bootstrapped.9p<-apply(bootstrapped.results.all.9p, 2, quantile, prob=c(0.025, 0.975))
#
# ## manage the data
# bootstropped.results.stat.9p <- t(rbind(t(mean.bootstrapped.9p), t(sd.bootstrapped.9p)))
# colnames(bootstropped.results.stat.9p) <- c("Mean", "sd")
# bootstropped.results.stat.9p <- cbind(bootstropped.results.stat.9p, t(qt.bootstrapped.9p))
#
# write.csv(bootstropped.results.stat.9p, file = "Bootstrapped_statistic_9species_12Oct2017.csv")
#
# ##########################
# ### Pair-wise t test  ####
# ##########################
# ##### two habitats
# habitat.results.test.9p <- data.frame()
# habitat.results.test.9p<-t.test(bootstrapped.results.all.9p[,1], bootstrapped.results.all.9p[,2])
# habitat.results.test.9p<-rbind(c(habitat.results.test.9p$statistic,
#                               habitat.results.test.9p$parameter,
#                               habitat.results.test.9p$p.value))
# colnames(habitat.results.test.9p)<-c("t","df","P value")
# rownames(habitat.results.test.9p) <- habitat.pair
#
# ##### four climate
# climate.results.test.9p<-data.frame()
#
# for (z in 3:5)
# {
#   for (h in z:5) {
#     temp.result.test.9p<-t.test(bootstrapped.results.all.9p[,z], bootstrapped.results.all.9p[, h+1])
#     climate.results.test.9p<-rbind(climate.results.test.9p,c(temp.result.test.9p$statistic,temp.result.test.9p$parameter,temp.result.test.9p$p.value))
#   }
# }
#
# colnames(climate.results.test.9p)<-c("t","df","P value")
# rownames(climate.results.test.9p) <- climate.pair
#
# #### eight vegetations
# vegetation.results.test.9p<-data.frame()
#
# for (w in 7:13)
# {
#   for (q in w:13) {
#     temp.result.test.9p<-t.test(bootstrapped.results.all.9p[,w], bootstrapped.results.all.9p[, q+1])
#     vegetation.results.test.9p<-rbind(vegetation.results.test.9p,c(temp.result.test.9p$statistic,temp.result.test.9p$parameter,temp.result.test.9p$p.value))
#   }
# }
#
# colnames(vegetation.results.test.9p)<-c("t","df","P value")
# rownames(vegetation.results.test.9p) <- vegetation.pair
#
# ## combine all results
# t.test.results.all.9species <- rbind(habitat.results.test.9p,
#                                      climate.results.test.9p,
#                                      vegetation.results.test.9p)
#
# write.csv(t.test.results.all.9species, file = "Bootstrapped_statistic_table_9species_11Oct2017.csv")
# #
# # ##### t.test with TJS (7-species community)
# # t.test.results.TJS<-data.frame()
# #
# # for (i in 1:14) {
# #   t.test.results.TJS.temp<-t.test(bootstrapped.results.all[,i], alternative="two.sided", mu = mean(eds.tjs[,1]))
# #   t.test.results.TJS<-rbind(t.test.results.TJS, c(t.test.results.TJS.temp$statistic,
# #                                                     t.test.results.TJS.temp$parameter,
# #                                                     t.test.results.TJS.temp$p.value))
# # }
# #
# # colnames(t.test.results.TJS)<-c("t","df","P value")
# # # write.csv(t.test.results.TJS, file = "Bootstrapped statistic table TJS 7_species.csv")
#
# ## histgram plots
# par(mfcol=c(8,3),
#     oma=c(1,1,0,1),
#     mar=c(2,2,1,0),
#     mgp=c(1,0.6,0),
#     tck=-0.02)
#
# # 2 habitats
# hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
# hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
# hist(cl.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
#
# h<-hist(cl.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey30")
# xfit<-seq(min(cl.results[,1]), max(cl.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(cl.results[,1]),sd=sd(cl.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(cl.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(cl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(op.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,120), ann=F, col="black")
# xfit<-seq(min(op.results[,1]), max(op.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(op.results[,1]),sd=sd(op.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(op.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(op.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
# hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
# hist(op.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120), ann=F) ## occupy the plot space
#
# # 4 climates
# hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
# hist(tc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
#
# h<-hist(tc.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,130), ann=F, col="darkorange")
# xfit<-seq(min(tc.results[,1]), max(tc.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(tc.results[,1]),sd=sd(tc.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(tc.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(tc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
#
# h<-hist(ac.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="grey20")
# xfit<-seq(min(ac.results[,1]), max(ac.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(ac.results[,1]),sd=sd(ac.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(ac.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(ac.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(ttc.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
# xfit<-seq(min(ttc.results[,1]), max(ttc.results[,1]), length=40)
# yfit<-dnorm(xfit, mean=mean(ttc.results[,1]),sd=sd(ttc.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(ttc.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(ttc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(cc.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,120),
#         xlab="Bootstrapped ecological disparity", ylab="", main="", col="dodgerblue3")
# xfit<-seq(min(cc.results[,1]), max(cc.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(cc.results[,1]),sd=sd(cc.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(cc.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(cc.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
# hist(cc.results[,1], breaks = 40, xlim = c(0,8), ylim = c(0,120),ann=F) ## occupy the plot space
#
# ## 8 vegetations
# h<-hist(trf.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,140), ann=F, col="red")
# xfit<-seq(min(trf.results[,1]), max(trf.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(trf.results[,1]),sd=sd(trf.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(trf.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(trf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(tsf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="orangered")
# xfit<-seq(min(tsf.results[,1]), max(tsf.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(tsf.results[,1]),sd=sd(tsf.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(tsf.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(tsf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(sv.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="firebrick1")
# xfit<-seq(min(sv.results[,1]), max(sv.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(sv.results[,1]),sd=sd(sv.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(sv.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(sv.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(gl.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,140), ann=F, col="goldenrod")
# xfit<-seq(min(gl.results[,1]), max(gl.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(gl.results[,1]),sd=sd(gl.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(gl.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(gl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(sl.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkgoldenrod1")
# xfit<-seq(min(sl.results[,1]), max(sl.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(sl.results[,1]),sd=sd(sl.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(sl.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(sl.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(dt.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,200), ann=F, col= "khaki")
# xfit<-seq(min(dt.results[,1]), max(dt.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(dt.results[,1]),sd=sd(dt.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(dt.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(dt.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(ttf.results[,1], breaks = 30, xlim = c(0,8), ylim = c(0,120), ann=F, col="green3")
# xfit<-seq(min(ttf.results[,1]), max(ttf.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(ttf.results[,1]),sd=sd(ttf.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(ttf.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(ttf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
# h<-hist(bf.results[,1], breaks = 20, xlim = c(0,8), ylim = c(0,120), ann=F, col="darkviolet")
# xfit<-seq(min(bf.results[,1]), max(bf.results[,1]), length=40)
# yfit<-dnorm(xfit,mean=mean(bf.results[,1]),sd=sd(bf.results[,1]))
# yfit <- yfit*diff(h$mids[1:2])*length(bf.results[,1])
# lines(xfit, yfit, col="blue", lwd=1)
# abline(v=quantile(bf.results[,1], prob=c(0.025, 0.975)), col="red", lwd=1.5, lty=2)
#
