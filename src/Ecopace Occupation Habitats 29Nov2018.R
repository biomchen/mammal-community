### The script code here is used to test whether the alternative placement of the saltotorial locomotor mode would drastically
### change the implications of the results

### Possible libraries used in the data analyses
library(permute)
library(lattice)
library(vegan)
library(MASS)
library(lattice)
library(boot)
library(sp)
library(parallel)
library(gridExtra)
library(ggplot2)

setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis")

### Read the database of the 100 communities
mcm.habitat.alt<-read.csv(file="EcoHabitatAll_6June2017.csv", sep=",", header=T)

## identify the russia park, four chinese sites 68-72, wetlands, and overgrassed, as well communities 102 and 105
rownames(mcm.habitat.alt[(mcm.habitat.alt$Community.No == 68),])
rownames(mcm.habitat.alt[(mcm.habitat.alt$Community.No == 72),])
rownames(mcm.habitat.alt[(mcm.habitat.alt$Community.No == 102),])
rownames(mcm.habitat.alt[(mcm.habitat.alt$Community.No == 105),])

## remove those unwanted communities
mcm.habitat.alt<- mcm.habitat.alt[-c(546:573, 821:823, 830:831),]

## save the current data to a csv file
write.csv(mcm.habitat.alt, file="EcoHabitatAll_Alt_29Nov2018.csv")

### the number of species in each commmunity
n.species.each.community.alt <- table(mcm.habitat.alt$Community.No)
n.species.each.community.alt <- as.data.frame(n.species.each.community.alt)
colnames(n.species.each.community.alt) <- c("Community", "No_species")

### alternative chocie of the locomotor ranking
### gliding=1, arboreal=2, scansorial=3, terrestrial=4, saltatorial=5, semiaquaitc=6, semifossorial=7, fossorial=8
### make it as string for replace function
str(mcm.habitat.alt)

# change their ranks
# scheme a
# # replace saltatorial=8 to saltorail=9
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==8,9)
# # replace semiaquatic=5 to 6
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==5,6)
# # replace semifossorial=6 to 7
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==6,7)
# # replace fossorial=7 to 8
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==7,8)
# # change back saltatorial=9 to saltorail=5
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==9,5)
# 
# # scheme b
# # replace saltatorial=8 to saltorail=9
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==8,9)
# # replace semifossorial=6 to 7
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==6,7)
# # replace fossorial=7 to 8
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==7,8)
# # change back saltatorial=9 to saltorail=5
# mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==9,6)

# scheme c
# replace saltatorial=8 to saltorail=9
mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==8,9)
# replace fossorial=7 to 8
mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==7,8)
# change back saltatorial=9 to saltorail=5
mcm.habitat.alt$Locomotor.Rank<-replace (mcm.habitat.alt$Locomotor.Rank, mcm.habitat.alt$Locomotor.Rank==9,7)

#########################################################
##### EDisp + ERich Aanalyses based on two different habitat types ####
#########################################################

#### open habitat
mcm.open<-mcm.habitat.alt[mcm.habitat.alt$Habitat == "Open",]

#### close habitat
mcm.close<-mcm.habitat.alt[mcm.habitat.alt$Habitat == "Close",]

## open habitat
ed.open<-data.frame()
dep.open<-data.frame()
cn.open<-unique(mcm.open$Community.No) ## The number of commuities
ed.open.data<-data.frame()

## disparity parameter of body size
bs.temp.open<-data.frame()
bs.community.open<-data.frame()
bs.open<-data.frame()
mean.bs.community.open<-data.frame()

## disparity parameter of diet
dt.temp.open<-data.frame()
dt.community.open<-data.frame()
dt.open<-data.frame()
mean.dt.community.open<-data.frame()

## disparity parameter of locomotor mode
lm.temp.open<-data.frame()
lm.community.open<-data.frame()
lm.open<-data.frame()
mean.lm.community.open<-data.frame()


for (j.open in cn.open) {
  newdata.open<-mcm.open[mcm.open$Community.No==j.open,]
  ed.open.data<-rbind(ed.open.data,c(newdata.open[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.open<-newdata.open[c(3,4,5)]
  result.open<-nrow(unique(data.open))
  ed.open<-rbind(ed.open,c(result.open))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.open in 1:(dim(data.open)[1]-1)) {
    for (v.open in u.open:(dim(data.open)[1]-1)) {
      ### disparity
      results.temp.open<-abs(data.open[u.open,]-data.open[v.open+1,])
      dep.temp.open<-rowSums(results.temp.open)
      dep.open<-rbind(dep.open,c(dep.temp.open))
    }
  }
}


## Close habitat
ed.close<-data.frame()
dep.close<-data.frame()
cn.close<-unique(mcm.close$Community.No) ## The number of commuities
ed.close.data<-data.frame()

## disparity parameter of body size
bs.temp.close<-data.frame()
bs.community.close<-data.frame()
bs.close<-data.frame()
mean.bs.community.close<-data.frame()

## disparity parameter of diet
dt.temp.close<-data.frame()
dt.community.close<-data.frame()
dt.close<-data.frame()
mean.dt.community.close<-data.frame()

## disparity parameter of locomotor mode
lm.temp.close<-data.frame()
lm.community.close<-data.frame()
lm.close<-data.frame()
mean.lm.community.close<-data.frame()


for (j.close in cn.close) {
  newdata.close<-mcm.close[mcm.close$Community.No==j.close,]
  ed.close.data<-rbind(ed.close.data,c(newdata.close[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.close<-newdata.close[c(3,4,5)]
  result.close<-nrow(unique(data.close))
  ed.close<-rbind(ed.close,c(result.close))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.close in 1:(dim(data.close)[1]-1)) {
    for (v.close in u.close:(dim(data.close)[1]-1)) {
      ### disparity
      results.temp.close<-abs(data.close[u.close,]-data.close[v.close+1,])
      dep.temp.close<-rowSums(results.temp.close)
      dep.close<-rbind(dep.close,c(dep.temp.close))
    }
  }
}

### EDisp
habitat.dep.mean <- rbind(mean(dep.close[,1]),
                          mean(dep.open[,1]))
habitat.dep.sd <- rbind(sd(dep.close[,1]),
                        sd(dep.open[,1]))
habitat.dep <- cbind(habitat.dep.mean, habitat.dep.sd)
colnames(habitat.dep) <- c("Mean", "SD")
rownames(habitat.dep) <- c("Closed", "Open")

### ERich
habitat.ed.mean <- rbind(mean(ed.close[,1]),
                         mean(ed.open[,1]))
habitat.ed.sd <- rbind(sd(ed.close[,1]),
                       sd(ed.open[,1]))
habitat.ed <- cbind(habitat.ed.mean, habitat.ed.sd)
colnames(habitat.ed) <- c("Mean", "SD")
rownames(habitat.ed) <- c("Closed", "Open")


###############################################################
##### EDisp + Erich based on four different climatic types ####
###############################################################

mcm.tropical.climate<-mcm.habitat.alt[mcm.habitat.alt$Climate == "Tropical",]

mcm.arid.climate<-mcm.habitat.alt[mcm.habitat.alt$Climate == "Arid",]

mcm.temperate.climate<-mcm.habitat.alt[mcm.habitat.alt$Climate== "Temperate",]

mcm.cold.climate<-mcm.habitat.alt[mcm.habitat.alt$Climate == "Cold",]

## tropical region
ed.tropical.climate<-data.frame()
dep.tropical.climate<-data.frame()
cn.tropical.climate<-unique(mcm.tropical.climate$Community.No) ## The number of commuities
ed.tropical.climate.data<-data.frame()

## disparity parameter of body size
bs.temp.tropical.climate<-data.frame()
bs.community.tropical.climate<-data.frame()
bs.tropical.climate<-data.frame()
mean.bs.community.tropical.climate<-data.frame()

## disparity parameter of diet
dt.temp.tropical.climate<-data.frame()
dt.community.tropical.climate<-data.frame()
dt.tropical.climate<-data.frame()
mean.dt.community.tropical.climate<-data.frame()

## disparity parameter of locomotor mode
lm.temp.tropical.climate<-data.frame()
lm.community.tropical.climate<-data.frame()
lm.tropical.climate<-data.frame()
mean.lm.community.tropical.climate<-data.frame()


for (j.tropical.climate in cn.tropical.climate) {
  newdata.tropical.climate<-mcm.tropical.climate[mcm.tropical.climate$Community.No==j.tropical.climate,]
  ed.tropical.climate.data<-rbind(ed.tropical.climate.data,c(newdata.tropical.climate[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.tropical.climate<-newdata.tropical.climate[c(3,4,5)]
  result.tropical.climate<-nrow(unique(data.tropical.climate))
  ed.tropical.climate<-rbind(ed.tropical.climate,c(result.tropical.climate))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.tropical.climate in 1:(dim(data.tropical.climate)[1]-1)) {
    for (v.tropical.climate in u.tropical.climate:(dim(data.tropical.climate)[1]-1)) {
      ### disparity
      results.temp.tropical.climate<-abs(data.tropical.climate[u.tropical.climate,]-data.tropical.climate[v.tropical.climate+1,])
      dep.temp.tropical.climate<-rowSums(results.temp.tropical.climate)
      dep.tropical.climate<-rbind(dep.tropical.climate,c(dep.temp.tropical.climate))
    }
  }
}


## arid region
ed.arid.climate<-data.frame()
dep.arid.climate<-data.frame()
cn.arid.climate<-unique(mcm.arid.climate$Community.No) ## The number of commuities
ed.arid.climate.data<-data.frame()

## disparity parameter of body size
bs.temp.arid.climate<-data.frame()
bs.community.arid.climate<-data.frame()
bs.arid.climate<-data.frame()
mean.bs.community.arid.climate<-data.frame()

## disparity parameter of diet
dt.temp.arid.climate<-data.frame()
dt.community.arid.climate<-data.frame()
dt.arid.climate<-data.frame()
mean.dt.community.arid.climate<-data.frame()

## disparity parameter of locomotor mode
lm.temp.arid.climate<-data.frame()
lm.community.arid.climate<-data.frame()
lm.arid.climate<-data.frame()
mean.lm.community.arid.climate<-data.frame()

for (j.arid.climate in cn.arid.climate) {
  newdata.arid.climate<-mcm.arid.climate[mcm.arid.climate$Community.No==j.arid.climate,]
  ed.arid.climate.data<-rbind(ed.arid.climate.data,c(newdata.arid.climate[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.arid.climate<-newdata.arid.climate[c(3,4,5)]
  result.arid.climate<-nrow(unique(data.arid.climate))
  ed.arid.climate<-rbind(ed.arid.climate,c(result.arid.climate))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.arid.climate in 1:(dim(data.arid.climate)[1]-1)) {
    for (v.arid.climate in u.arid.climate:(dim(data.arid.climate)[1]-1)) {
      ### disparity
      results.temp.arid.climate<-abs(data.arid.climate[u.arid.climate,]-data.arid.climate[v.arid.climate+1,])
      dep.temp.arid.climate<-rowSums(results.temp.arid.climate)
      dep.arid.climate<-rbind(dep.arid.climate,c(dep.temp.arid.climate))
    }
  }
}


## temperate region
ed.temperate.climate<-data.frame()
dep.temperate.climate<-data.frame()
cn.temperate.climate<-unique(mcm.temperate.climate$Community.No) ## The number of commuities
ed.temperate.climate.data<-data.frame()

## disparity parameter of body size
bs.temp.temperate.climate<-data.frame()
bs.community.temperate.climate<-data.frame()
bs.temperate.climate<-data.frame()
mean.bs.community.temperate.climate<-data.frame()

## disparity parameter of diet
dt.temp.temperate.climate<-data.frame()
dt.community.temperate.climate<-data.frame()
dt.temperate.climate<-data.frame()
mean.dt.community.temperate.climate<-data.frame()

## disparity parameter of locomotor mode
lm.temp.temperate.climate<-data.frame()
lm.community.temperate.climate<-data.frame()
lm.temperate.climate<-data.frame()
mean.lm.community.temperate.climate<-data.frame()


for (j.temperate.climate in cn.temperate.climate) {
  newdata.temperate.climate<-mcm.temperate.climate[mcm.temperate.climate$Community.No==j.temperate.climate,]
  ed.temperate.climate.data<-rbind(ed.temperate.climate.data,c(newdata.temperate.climate[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.temperate.climate<-newdata.temperate.climate[c(3,4,5)]
  result.temperate.climate<-nrow(unique(data.temperate.climate))
  ed.temperate.climate<-rbind(ed.temperate.climate,c(result.temperate.climate))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.temperate.climate in 1:(dim(data.temperate.climate)[1]-1)) {
    for (v.temperate.climate in u.temperate.climate:(dim(data.temperate.climate)[1]-1)) {
      ### disparity
      results.temp.temperate.climate<-abs(data.temperate.climate[u.temperate.climate,]-data.temperate.climate[v.temperate.climate+1,])
      dep.temp.temperate.climate<-rowSums(results.temp.temperate.climate)
      dep.temperate.climate<-rbind(dep.temperate.climate,c(dep.temp.temperate.climate))
    }
  }
}

## cold region
ed.cold.climate<-data.frame()
dep.cold.climate<-data.frame()
cn.cold.climate<-unique(mcm.cold.climate$Community.No) ## The number of commuities
ed.cold.climate.data<-data.frame()

## disparity parameter of body size
bs.temp.cold.climate<-data.frame()
bs.community.cold.climate<-data.frame()
bs.cold.climate<-data.frame()
mean.bs.community.cold.climate<-data.frame()

## disparity parameter of diet
dt.temp.cold.climate<-data.frame()
dt.community.cold.climate<-data.frame()
dt.cold.climate<-data.frame()
mean.dt.community.cold.climate<-data.frame()

## disparity parameter of locomotor mode
lm.temp.cold.climate<-data.frame()
lm.community.cold.climate<-data.frame()
lm.cold.climate<-data.frame()
mean.lm.community.cold.climate<-data.frame()


for (j.cold.climate in cn.cold.climate) {
  newdata.cold.climate<-mcm.cold.climate[mcm.cold.climate$Community.No==j.cold.climate,]
  ed.cold.climate.data<-rbind(ed.cold.climate.data,c(newdata.cold.climate[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.cold.climate<-newdata.cold.climate[c(3,4,5)]
  result.cold.climate<-nrow(unique(data.cold.climate))
  ed.cold.climate<-rbind(ed.cold.climate,c(result.cold.climate))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.cold.climate in 1:(dim(data.cold.climate)[1]-1)) {
    for (v.cold.climate in u.cold.climate:(dim(data.cold.climate)[1]-1)) {
      ### disparity
      results.temp.cold.climate<-abs(data.cold.climate[u.cold.climate,]-data.cold.climate[v.cold.climate+1,])
      dep.temp.cold.climate<-rowSums(results.temp.cold.climate)
      dep.cold.climate<-rbind(dep.cold.climate,c(dep.temp.cold.climate))
    }
  }
}


### EDisp
climate.dep.mean <- rbind(mean(dep.tropical.climate[,1]),
                          mean(dep.arid.climate[,1]),
                          mean(dep.temperate.climate[,1]),
                          mean(dep.cold.climate[,1]))
climate.dep.sd <- rbind(sd(dep.tropical.climate[,1]),
                        sd(dep.arid.climate[,1]),
                        sd(dep.temperate.climate[,1]),
                        sd(dep.cold.climate[,1]))
climate.dep <- cbind(climate.dep.mean, climate.dep.sd)
colnames(climate.dep) <- c("Mean", "SD")
rownames(climate.dep) <- c("Tropical", "Arid","Temperate", "Cold")

### ERich
climate.ed.mean <- rbind(mean(ed.tropical.climate[,1]),
                         mean(ed.arid.climate[,1]),
                         mean(ed.temperate.climate[,1]),
                         mean(ed.cold.climate[,1]))
climate.ed.sd <- rbind(sd(ed.tropical.climate[,1]),
                       sd(ed.arid.climate[,1]),
                       sd(ed.temperate.climate[,1]),
                       sd(ed.cold.climate[,1]))
climate.ed <- cbind(climate.ed.mean, climate.ed.sd)
colnames(climate.ed) <- c("Mean", "SD")
rownames(climate.ed) <- c("Tropical", "Arid","Temperate", "Cold")


####################################################
##### EDisp + Erich based on eight vegetation types ####
####################################################

##### The data of each habitat
# Boreal forest
mcm.boreal.forest<-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Boreal forest",]
# Desert
mcm.desert <-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Desert",]
# Temperate forest
mcm.temperate.forest<-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Temperate forest",]
# Grassland
mcm.grassland<-mcm.habitat.alt[mcm.habitat.alt$Vegetation== "Grassland",]
# Shrubland
mcm.shrubland<-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Shrubland",]
# Savanna
mcm.savanna<-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Savanna",]
# Tropical seasonal forest
mcm.tropical.seasonal.forest<-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Tropical seasonal forest",]
# Tropical rainforest
mcm.tropical.rainforest<-mcm.habitat.alt[mcm.habitat.alt$Vegetation == "Tropical rainforest",]

### Tropical rainforest
ed.tropical.rainforest<-data.frame()
dep.tropical.rainforest<-data.frame()
cn.tropical.rainforest<-unique(mcm.tropical.rainforest$Community.No) ## The number of commuities
ed.tropical.rainforest.data<-data.frame()

## disparity parameter of body size
bs.temp.tropical.rainforest<-data.frame()
bs.community.tropical.rainforest<-data.frame()
bs.tropical.rainforest<-data.frame()
mean.bs.community.tropical.rainforest<-data.frame()

## disparity parameter of diet
dt.temp.tropical.rainforest<-data.frame()
dt.community.tropical.rainforest<-data.frame()
dt.tropical.rainforest<-data.frame()
mean.dt.community.tropical.rainforest<-data.frame()

## disparity parameter of locomotor mode
lm.temp.tropical.rainforest<-data.frame()
lm.community.tropical.rainforest<-data.frame()
lm.tropical.rainforest<-data.frame()
mean.lm.community.tropical.rainforest<-data.frame()


for (j.tropical.rainforest in cn.tropical.rainforest) {
  newdata.tropical.rainforest<-mcm.tropical.rainforest[mcm.tropical.rainforest$Community.No==j.tropical.rainforest,]
  ed.tropical.rainforest.data<-rbind(ed.tropical.rainforest.data,c(newdata.tropical.rainforest[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.tropical.rainforest<-newdata.tropical.rainforest[c(3,4,5)]
  result.tropical.rainforest<-nrow(unique(data.tropical.rainforest))
  ed.tropical.rainforest<-rbind(ed.tropical.rainforest,c(result.tropical.rainforest))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.tropical.rainforest in 1:(dim(data.tropical.rainforest)[1]-1)) {
    for (v.tropical.rainforest in u.tropical.rainforest:(dim(data.tropical.rainforest)[1]-1)) {
      ### disparity
      results.temp.tropical.rainforest<-abs(data.tropical.rainforest[u.tropical.rainforest,]-data.tropical.rainforest[v.tropical.rainforest+1,])
      dep.temp.tropical.rainforest<-rowSums(results.temp.tropical.rainforest)
      dep.tropical.rainforest<-rbind(dep.tropical.rainforest,c(dep.temp.tropical.rainforest))
    }
  }
}



### Tropical seasonal forest
ed.tropical.seasonal.forest<-data.frame()
dep.tropical.seasonal.forest<-data.frame()
cn.tropical.seasonal.forest<-unique(mcm.tropical.seasonal.forest$Community.No) ## The number of commuities
ed.tropical.seasonal.forest.data<-data.frame()

## disparity parameter of body size
bs.temp.tropical.seasonal.forest<-data.frame()
bs.community.tropical.seasonal.forest<-data.frame()
bs.tropical.seasonal.forest<-data.frame()
mean.bs.community.tropical.seasonal.forest<-data.frame()

## disparity parameter of diet
dt.temp.tropical.seasonal.forest<-data.frame()
dt.community.tropical.seasonal.forest<-data.frame()
dt.tropical.seasonal.forest<-data.frame()
mean.dt.community.tropical.seasonal.forest<-data.frame()

## disparity parameter of locomotor mode
lm.temp.tropical.seasonal.forest<-data.frame()
lm.community.tropical.seasonal.forest<-data.frame()
lm.tropical.seasonal.forest<-data.frame()
mean.lm.community.tropical.seasonal.forest<-data.frame()



for (j.tropical.seasonal.forest in cn.tropical.seasonal.forest) {
  newdata.tropical.seasonal.forest<-mcm.tropical.seasonal.forest[mcm.tropical.seasonal.forest$Community.No==j.tropical.seasonal.forest,]
  ed.tropical.seasonal.forest.data<-rbind(ed.tropical.seasonal.forest.data,c(newdata.tropical.seasonal.forest[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.tropical.seasonal.forest<-newdata.tropical.seasonal.forest[c(3,4,5)]
  result.tropical.seasonal.forest<-nrow(unique(data.tropical.seasonal.forest))
  ed.tropical.seasonal.forest<-rbind(ed.tropical.seasonal.forest,c(result.tropical.seasonal.forest))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.tropical.seasonal.forest in 1:(dim(data.tropical.seasonal.forest)[1]-1)) {
    for (v.tropical.seasonal.forest in u.tropical.seasonal.forest:(dim(data.tropical.seasonal.forest)[1]-1)) {
      ### disparity
      results.temp.tropical.seasonal.forest<-abs(data.tropical.seasonal.forest[u.tropical.seasonal.forest,]-data.tropical.seasonal.forest[v.tropical.seasonal.forest+1,])
      dep.temp.tropical.seasonal.forest<-rowSums(results.temp.tropical.seasonal.forest)
      dep.tropical.seasonal.forest<-rbind(dep.tropical.seasonal.forest,c(dep.temp.tropical.seasonal.forest))
    }
  }
}

### savanna
ed.savanna<-data.frame()
dep.savanna<-data.frame()
cn.savanna<-unique(mcm.savanna$Community.No) ## The number of commuities
ed.savanna.data<-data.frame()

## disparity parameter of body size
bs.temp.savanna<-data.frame()
bs.community.savanna<-data.frame()
bs.savanna<-data.frame()
mean.bs.community.savanna<-data.frame()

## disparity parameter of diet
dt.temp.savanna<-data.frame()
dt.community.savanna<-data.frame()
dt.savanna<-data.frame()
mean.dt.community.savanna<-data.frame()

## disparity parameter of locomotor mode
lm.temp.savanna<-data.frame()
lm.community.savanna<-data.frame()
lm.savanna<-data.frame()
mean.lm.community.savanna<-data.frame()


for (j.savanna in cn.savanna) {
  newdata.savanna<-mcm.savanna[mcm.savanna$Community.No==j.savanna,]
  ed.savanna.data<-rbind(ed.savanna.data,c(newdata.savanna[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.savanna<-newdata.savanna[c(3,4,5)]
  result.savanna<-nrow(unique(data.savanna))
  ed.savanna<-rbind(ed.savanna,c(result.savanna))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.savanna in 1:(dim(data.savanna)[1]-1)) {
    for (v.savanna in u.savanna:(dim(data.savanna)[1]-1)) {
      ### disparity
      results.temp.savanna<-abs(data.savanna[u.savanna,]-data.savanna[v.savanna+1,])
      dep.temp.savanna<-rowSums(results.temp.savanna)
      dep.savanna<-rbind(dep.savanna,c(dep.temp.savanna))
    }
  }
}



### Grassland
ed.grassland<-data.frame()
dep.grassland<-data.frame()
cn.grassland<-unique(mcm.grassland$Community.No) ## The number of commuities
ed.grassland.data<-data.frame()

## disparity parameter of body size
bs.temp.grassland<-data.frame()
bs.community.grassland<-data.frame()
bs.grassland<-data.frame()
mean.bs.community.grassland<-data.frame()

## disparity parameter of diet
dt.temp.grassland<-data.frame()
dt.community.grassland<-data.frame()
dt.grassland<-data.frame()
mean.dt.community.grassland<-data.frame()

## disparity parameter of locomotor mode
lm.temp.grassland<-data.frame()
lm.community.grassland<-data.frame()
lm.grassland<-data.frame()
mean.lm.community.grassland<-data.frame()


for (j.grassland in cn.grassland) {
  newdata.grassland<-mcm.grassland[mcm.grassland$Community.No==j.grassland,]
  ed.grassland.data<-rbind(ed.grassland.data,c(newdata.grassland[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.grassland<-newdata.grassland[c(3,4,5)]
  result.grassland<-nrow(unique(data.grassland))
  ed.grassland<-rbind(ed.grassland,c(result.grassland))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.grassland in 1:(dim(data.grassland)[1]-1)) {
    for (v.grassland in u.grassland:(dim(data.grassland)[1]-1)) {
      ### disparity
      results.temp.grassland<-abs(data.grassland[u.grassland,]-data.grassland[v.grassland+1,])
      dep.temp.grassland<-rowSums(results.temp.grassland)
      dep.grassland<-rbind(dep.grassland,c(dep.temp.grassland))
    }
  }
}


### Shrubland
ed.shrubland<-data.frame()
dep.shrubland<-data.frame()
cn.shrubland<-unique(mcm.shrubland$Community.No) ## The number of commuities
ed.shrubland.data<-data.frame()

## disparity parameter of body size
bs.temp.shrubland<-data.frame()
bs.community.shrubland<-data.frame()
bs.shrubland<-data.frame()
mean.bs.community.shrubland<-data.frame()

## disparity parameter of diet
dt.temp.shrubland<-data.frame()
dt.community.shrubland<-data.frame()
dt.shrubland<-data.frame()
mean.dt.community.shrubland<-data.frame()

## disparity parameter of locomotor mode
lm.temp.shrubland<-data.frame()
lm.community.shrubland<-data.frame()
lm.shrubland<-data.frame()
mean.lm.community.shrubland<-data.frame()


for (j.shrubland in cn.shrubland) {
  newdata.shrubland<-mcm.shrubland[mcm.shrubland$Community.No==j.shrubland,]
  ed.shrubland.data<-rbind(ed.shrubland.data,c(newdata.shrubland[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.shrubland<-newdata.shrubland[c(3,4,5)]
  result.shrubland<-nrow(unique(data.shrubland))
  ed.shrubland<-rbind(ed.shrubland,c(result.shrubland))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.shrubland in 1:(dim(data.shrubland)[1]-1)) {
    for (v.shrubland in u.shrubland:(dim(data.shrubland)[1]-1)) {
      ### disparity
      results.temp.shrubland<-abs(data.shrubland[u.shrubland,]-data.shrubland[v.shrubland+1,])
      dep.temp.shrubland<-rowSums(results.temp.shrubland)
      dep.shrubland<-rbind(dep.shrubland,c(dep.temp.shrubland))
    }
  }
}

### Desert
ed.desert<-data.frame()
dep.desert<-data.frame()
cn.desert<-unique(mcm.desert$Community.No) ## The number of commuities
ed.desert.data<-data.frame()

## disparity parameter of body size
bs.temp.desert<-data.frame()
bs.community.desert<-data.frame()
bs.desert<-data.frame()
mean.bs.community.desert<-data.frame()

## disparity parameter of diet
dt.temp.desert<-data.frame()
dt.community.desert<-data.frame()
dt.desert<-data.frame()
mean.dt.community.desert<-data.frame()

## disparity parameter of locomotor mode
lm.temp.desert<-data.frame()
lm.community.desert<-data.frame()
lm.desert<-data.frame()
mean.lm.community.desert<-data.frame()


for (j.desert in cn.desert) {
  newdata.desert<-mcm.desert[mcm.desert$Community.No==j.desert,]
  ed.desert.data<-rbind(ed.desert.data,c(newdata.desert[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.desert<-newdata.desert[c(3,4,5)]
  result.desert<-nrow(unique(data.desert))
  ed.desert<-rbind(ed.desert,c(result.desert))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.desert in 1:(dim(data.desert)[1]-1)) {
    for (v.desert in u.desert:(dim(data.desert)[1]-1)) {
      ### disparity
      results.temp.desert<-abs(data.desert[u.desert,]-data.desert[v.desert+1,])
      dep.temp.desert<-rowSums(results.temp.desert)
      dep.desert<-rbind(dep.desert,c(dep.temp.desert))
    }
  }
}



### temperate forest
ed.temperate.forest<-data.frame()
dep.temperate.forest<-data.frame()
cn.temperate.forest<-unique(mcm.temperate.forest$Community.No) ## The number of commuities
ed.temperate.forest.data<-data.frame()

## disparity parameter of body size
bs.temp.temperate.forest<-data.frame()
bs.community.temperate.forest<-data.frame()
bs.temperate.forest<-data.frame()
mean.bs.community.temperate.forest<-data.frame()

## disparity parameter of diet
dt.temp.temperate.forest<-data.frame()
dt.community.temperate.forest<-data.frame()
dt.temperate.forest<-data.frame()
mean.dt.community.temperate.forest<-data.frame()

## disparity parameter of locomotor mode
lm.temp.temperate.forest<-data.frame()
lm.community.temperate.forest<-data.frame()
lm.temperate.forest<-data.frame()
mean.lm.community.temperate.forest<-data.frame()


for (j.temperate.forest in cn.temperate.forest) {
  newdata.temperate.forest<-mcm.temperate.forest[mcm.temperate.forest$Community.No==j.temperate.forest,]
  ed.temperate.forest.data<-rbind(ed.temperate.forest.data,c(newdata.temperate.forest[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.temperate.forest<-newdata.temperate.forest[c(3,4,5)]
  result.temperate.forest<-nrow(unique(data.temperate.forest))
  ed.temperate.forest<-rbind(ed.temperate.forest,c(result.temperate.forest))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.temperate.forest in 1:(dim(data.temperate.forest)[1]-1)) {
    for (v.temperate.forest in u.temperate.forest:(dim(data.temperate.forest)[1]-1)) {
      ### disparity
      results.temp.temperate.forest<-abs(data.temperate.forest[u.temperate.forest,]-data.temperate.forest[v.temperate.forest+1,])
      dep.temp.temperate.forest<-rowSums(results.temp.temperate.forest)
      dep.temperate.forest<-rbind(dep.temperate.forest,c(dep.temp.temperate.forest))
    }
  }
}

### boreal foest;
ed.boreal<-data.frame()
dep.boreal<-data.frame()
cn.boreal<-unique(mcm.boreal.forest$Community.No) ## The number of commuities
ed.boreal.data<-data.frame()

## disparity parameter of body size
bs.temp.boreal<-data.frame()
bs.community.boreal<-data.frame()
bs.boreal<-data.frame()
mean.bs.community.boreal<-data.frame()

## disparity parameter of diet
dt.temp.boreal<-data.frame()
dt.community.boreal<-data.frame()
dt.boreal<-data.frame()
mean.dt.community.boreal<-data.frame()

## disparity parameter of locomotor mode
lm.temp.boreal<-data.frame()
lm.community.boreal<-data.frame()
lm.boreal<-data.frame()
mean.lm.community.boreal<-data.frame()

for (j.boreal in cn.boreal) {
  newdata.boreal<-mcm.boreal.forest[mcm.boreal.forest$Community.No==j.boreal,]
  ed.boreal.data<-rbind(ed.boreal.data,c(newdata.boreal[c(3,4,5)])) ### entire ecological raw dataset
  ### Ecologoical diversity
  data.boreal<-newdata.boreal[c(3,4,5)]
  result.boreal<-nrow(unique(data.boreal))
  ed.boreal<-rbind(ed.boreal,c(result.boreal))   ## store the data in sequence
  ### Disparity and body size differences
  for (u.boreal in 1:(dim(data.boreal)[1]-1)) {
    for (v.boreal in u.boreal:(dim(data.boreal)[1]-1)) {
      ### disparity
      results.temp.boreal<-abs(data.boreal[u.boreal,]-data.boreal[v.boreal+1,])
      dep.temp.boreal<-rowSums(results.temp.boreal)
      dep.boreal<-rbind(dep.boreal,c(dep.temp.boreal))
    }
  }
}


### save the mean and sd of EDisp and ERich
### EDisp
vegetation.dep.mean <- rbind(mean(dep.tropical.rainforest[,1]),
                             mean(dep.tropical.seasonal.forest[,1]),
                             mean(dep.savanna[,1]),
                             mean(dep.grassland[,1]),
                             mean(dep.shrubland[,1]),
                             mean(dep.desert[,1]),
                             mean(dep.temperate.forest[,1]),
                             mean(dep.boreal[,1]))
vegetation.dep.sd <- rbind(sd(dep.tropical.rainforest[,1]),
                           sd(dep.tropical.seasonal.forest[,1]),
                           sd(dep.savanna[,1]),
                           sd(dep.grassland[,1]),
                           sd(dep.shrubland[,1]),
                           sd(dep.desert[,1]),
                           sd(dep.temperate.forest[,1]),
                           sd(dep.boreal[,1]))
vegetation.dep <- cbind(vegetation.dep.mean, vegetation.dep.sd)
colnames(vegetation.dep) <- c("Mean", "SD")
rownames(vegetation.dep) <- c("Tropical rainforest", "Tropical seasoanl forest",
                              "Savanna", "Grassland","Shrubland", "Desert", 
                              "Temperate forest", "Boreal forest")

### ERich
vegetation.ed.mean <- rbind(mean(ed.tropical.rainforest[,1]),
                            mean(ed.tropical.seasonal.forest[,1]),
                            mean(ed.savanna[,1]),
                            mean(ed.grassland[,1]),
                            mean(ed.shrubland[,1]),
                            mean(ed.desert[,1]),
                            mean(ed.temperate.forest[,1]),
                            mean(ed.boreal[,1]))
vegetation.ed.sd <- rbind(sd(ed.tropical.rainforest[,1]),
                          sd(ed.tropical.seasonal.forest[,1]),
                          sd(ed.savanna[,1]),
                          sd(ed.grassland[,1]),
                          sd(ed.shrubland[,1]),
                          sd(ed.desert[,1]),
                          sd(ed.temperate.forest[,1]),
                          sd(ed.boreal[,1]))
vegetation.ed <- cbind(vegetation.ed.mean, vegetation.ed.sd)
colnames(vegetation.ed) <- c("Mean", "SD")
rownames(vegetation.ed) <- c("Tropical rainforest", "Tropical seasoanl forest",
                             "Savanna", "Grassland","Shrubland", "Desert", 
                             "Temperate forest", "Boreal forest")

## report the results
habitat.dep
climate.dep
vegetation.dep

dep.all.alt <- rbind(habitat.dep, climate.dep, vegetation.dep)
write.csv(dep.all.alt, file ="EDisp_all_alt__SchemeC_3Dec2018.csv")