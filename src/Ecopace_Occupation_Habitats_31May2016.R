### Ecological diversity analyses codes originally written by Meng Chen during the doctoral dissertation
### Reorder method is using climate region based on climate codes

### The purpose of the analyses is to aim to calculate how many different combination of three ecological 
### paramters in each mammal community to decipher the ecological differences among them

##############################################################################################################
### Based on new data added during Postdoc at Smithosinan, some of the codes has been revised for better   ###
### analyzing the data and represeneitng the visiualazaiton of the results                                 ###
### Meng Chen @ copyright                                              Date: 05Oct2015   Updated12Jan2016  ###
##############################################################################################################

####################################################################################################################
# Tests used in this study: 1) investigate the relationships among three ecological parameters using Kendal's test
#                           2) calculate the ecological dispairty and diversity, and with bootstrapping
#                           3) calculate the dissmiliarity index, and with bootstrapping
#                           4) use the GLM to model the ecoloigcal and disparity and diversity to investigate the
#                                 environmental influences on small-bodied mamaml communities
####################################################################################################################                           

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
mcm.habitat<-read.csv(file="EcoHabitatAll_6June2017.csv", sep=",", header=T)

## identify the russia park, four chinese sites 68-72, wetlands, and overgrassed, as well communities 102 and 105
rownames(mcm.habitat[(mcm.habitat$Community.No == 68),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 72),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 102),])
rownames(mcm.habitat[(mcm.habitat$Community.No == 105),])

## remove those unwanted communities
mcm.habitat<- mcm.habitat[-c(546:573, 821:823, 830:831),]

### the number of species in each commmunity
n.species.each.community <- table(mcm.habitat$Community.No)
n.species.each.community <- as.data.frame(n.species.each.community)
colnames(n.species.each.community) <- c("Community", "No_species")

#############################
### Ecotype identification
############################

## ecotype data
Eco.type.data <- as.data.frame(with(mcm.habitat, paste0(BS.Rank, Diet.Rank, Locomotor.Rank)))
colnames(Eco.type.data) <- c("EcoType")
Eco.type.data <-cbind(mcm.habitat$Community.No, mcm.habitat$Habitat, mcm.habitat$Climate,
                      mcm.habitat$Vegetation, Eco.type.data)
colnames(Eco.type.data) <- c("Community", "Habitat", "Climate", "Vegetation", "Ecotype")

write.csv(Eco.type.data, file = "Eco.type.data_Oct262017.csv")
ecotype.summary <- table(Eco.type.data$Ecotype)
write.csv(ecotype.summary, file = "Ecotype_summary_19Sept2018.csv")
unique.ecotype.global <- rownames(ecotype.summary)
unique.ecotype.global <- data.frame(unique.ecotype.global)
colnames(unique.ecotype.global) <- c("ecotype")
unique.ecotype.global.splitted <- t(sapply(unique.ecotype.global$ecotype, 
                                           function(x) substring (x, 1:3, 1:3)))
unique.ecotype.global.splitted <- data.frame(unique.ecotype.global.splitted)
write.csv(unique.ecotype.global.splitted, file = "Unique_ecotype_global_20Sept2018.csv")

community.occurrence.data <- table(Eco.type.data[,c(1,5)])
write.csv(community.occurrence.data, file = "Community_occurrence_19Sept2019.csv")

## Ecotype by different evnironments
format.function <- function (input) {
  input <- data.frame(input)
  input$Ecocell <- rownames(input)
  colnames(input) <- c("Frequency", "Eco-cell")
  return(input)
}

## habitat openness
## open habitat
Eco.type.open.data <- Eco.type.data[Eco.type.data$Habitat=="Open",]
summary.eco.type.open <- summary(Eco.type.open.data$Ecotype)
summary.eco.type.open <- format.function(summary.eco.type.open)
write.csv(summary.eco.type.open, file = "summary.eco.type.open.csv")
## closed habitat
Eco.type.closed.data <- Eco.type.data[Eco.type.data$Habitat=="Close",]
summary.eco.type.closed <- summary(Eco.type.closed.data$Ecotype)
summary.eco.type.closed <- format.function(summary.eco.type.closed)
write.csv(summary.eco.type.closed, file = "summary.eco.type.closed.csv")

## four climates
## tropical climate
Eco.type.tropical.data <- Eco.type.data[Eco.type.data$Climate=="Tropical",]
summary.eco.type.tropical <- summary(Eco.type.tropical.data$Ecotype)
summary.eco.type.tropical <- format.function(summary.eco.type.tropical)
write.csv(summary.eco.type.tropical, file = "summary.eco.type.tropical.csv")
## arid climate
Eco.type.arid.data <- Eco.type.data[Eco.type.data$Climate=="Arid",]
summary.eco.type.arid <- summary(Eco.type.arid.data$Ecotype)
summary.eco.type.arid <- format.function(summary.eco.type.arid)
write.csv(summary.eco.type.arid, file = "summary.eco.type.arid.csv")
## temperate climate
Eco.type.temperate.data <- Eco.type.data[Eco.type.data$Climate=="Temperate",]
summary.eco.type.temperate <- summary(Eco.type.temperate.data$Ecotype)
summary.eco.type.temperate <- format.function(summary.eco.type.temperate)
write.csv(summary.eco.type.temperate, file = "summary.eco.type.temperate.csv")
## cold climate
Eco.type.cold.data <- Eco.type.data[Eco.type.data$Climate=="Cold",]
summary.eco.type.cold<- summary(Eco.type.cold.data$Ecotype)
summary.eco.type.cold <- format.function(summary.eco.type.cold)
write.csv(summary.eco.type.cold, file = "summary.eco.type.cold.csv")

## eight vegetations
## tropical rainforest
Eco.type.tropicalrainforest.data <- Eco.type.data[Eco.type.data$Vegetation=="Tropical rainforest",]
summary.eco.type.tropicalrainforest<- summary(Eco.type.tropicalrainforest.data$Ecotype)
summary.eco.type.tropicalrainforest <- format.function(summary.eco.type.tropicalrainforest)
write.csv(summary.eco.type.tropicalrainforest, file = "summary.eco.type.tropicalrainforest.csv")
## tropical seasonal forest
Eco.type.tropicalseasonalforest.data <- Eco.type.data[Eco.type.data$Vegetation=="Tropical seasonal forest",]
summary.eco.type.tropicalseasonalforest<- summary(Eco.type.tropicalseasonalforest.data$Ecotype)
summary.eco.type.tropicalseasonalforest <- format.function(summary.eco.type.tropicalseasonalforest)
write.csv(summary.eco.type.tropicalseasonalforest, file = "summary.eco.type.tropicalseasonalforest.csv")
## savanna
Eco.type.savanna.data <- Eco.type.data[Eco.type.data$Vegetation=="Savanna",]
summary.eco.type.savanna<- summary(Eco.type.savanna.data$Ecotype)
summary.eco.type.savanna <- format.function(summary.eco.type.savanna)
write.csv(summary.eco.type.savanna, file = "summary.eco.type.savanna.csv")
## grassland
Eco.type.grassland.data <- Eco.type.data[Eco.type.data$Vegetation=="Grassland",]
summary.eco.type.grassland<- summary(Eco.type.grassland.data$Ecotype)
summary.eco.type.grassland <- format.function(summary.eco.type.grassland)
write.csv(summary.eco.type.grassland, file = "summary.eco.type.grassland.csv")
## shrubland
Eco.type.shrubland.data <- Eco.type.data[Eco.type.data$Vegetation=="Shrubland",]
summary.eco.type.shrubland<- summary(Eco.type.shrubland.data$Ecotype)
summary.eco.type.shrubland <- format.function(summary.eco.type.shrubland)
write.csv(summary.eco.type.shrubland, file = "summary.eco.type.shrubland.csv")
## desert
Eco.type.desert.data <- Eco.type.data[Eco.type.data$Vegetation=="Desert",]
summary.eco.type.desert<- summary(Eco.type.desert.data$Ecotype)
summary.eco.type.desert <- format.function(summary.eco.type.desert)
write.csv(summary.eco.type.desert, file = "summary.eco.type.desert.csv")
## temperate forest
Eco.type.temperateforest.data <- Eco.type.data[Eco.type.data$Vegetation=="Temperate forest",]
summary.eco.type.temperateforest<- summary(Eco.type.temperateforest.data$Ecotype)
summary.eco.type.temperateforest <- format.function(summary.eco.type.temperateforest)
write.csv(summary.eco.type.temperateforest, file = "summary.eco.type.temperateforest.csv")
## boreal forest
Eco.type.borealforest.data <- Eco.type.data[Eco.type.data$Vegetation=="Boreal forest",]
summary.eco.type.borealforest<- summary(Eco.type.borealforest.data$Ecotype)
summary.eco.type.borealforest <- format.function(summary.eco.type.borealforest)
write.csv(summary.eco.type.borealforest, file = "summary.eco.type.borealforest.csv")


##############################################
## Jaccard Dissimilarity test pair-wise global communtiies
##############################################

# Jaccard disimilarity function
jaccard.dis<-data.frame()
dissimilarity.results <-data.frame()

Jaccard.function<-function (input.data) {
  # calculate the dissimilarity
  Community<-unique(input.data$Community)

  for (g in 1:(length(Community)-1)) {
    for (h in g:(length(Community)-1)) {
      func.community.temp1<-input.data[input.data$Community==Community[g],]
      func.community.temp2<-unique(func.community.temp1$Ecotype)
      func.community.temp3<-input.data[input.data$Community==Community[h+1],]
      func.community.temp4<-unique(func.community.temp3$Ecotype)
      # jaccarad disimiliarity equation
      jaccard.dis.temp<-1-length(intersect(func.community.temp2,func.community.temp4))/
        length(union(func.community.temp2,func.community.temp4))
      jaccard.dis<-rbind(jaccard.dis,c(jaccard.dis.temp))
    }
  }
  return(jaccard.dis)
}

dissimilarity.results <- Jaccard.function(Eco.type.data)


########################################################
### proportion transoframtion of the occurrence data
########################################################

### ecoloigcal parameter occurrence data
mcm.habitat.occurrency<-cbind(table(mcm.habitat[,c(1,3)]), 
                              table(mcm.habitat[,c(1,4)]), 
                              table(mcm.habitat[,c(1,5)]))

colnames(mcm.habitat.occurrency) <- c("BS1", "BS2", "BS3", "BS4", "BS5",
                                      "DP1", "DP2", "DP3", "DP4", "DP5", "DP6",
                                      "LM1", "LM2", "LM3", "LM4", "LM5", "LM6", "LM7", "LM8")

#write.csv(mcm.habitat.occurrency, file = "mcm.occurrency.original_29May2017.csv")

### proportion transoframtion
mcm.habitat.proportion<-cbind(table(mcm.habitat[,c(1,3)])/rowSums(table(mcm.habitat[,c(1,3)])), 
                              table(mcm.habitat[,c(1,4)])/rowSums(table(mcm.habitat[,c(1,4)])), 
                              table(mcm.habitat[,c(1,5)])/rowSums(table(mcm.habitat[,c(1,5)])))

colnames(mcm.habitat.proportion) <- c("BS1", "BS2", "BS3", "BS4", "BS5",
                                      "DP1", "DP2", "DP3", "DP4", "DP5", "DP6",
                                      "LM1", "LM2", "LM3", "LM4", "LM5", "LM6", "LM7", "LM8")

#write.csv(mcm.habitat.proportion, file = "mcm.proportion.original_29May2017.csv")

############################################################
## investigate the cliamte and habitat independent variables
############################################################

tbl<-table(mcm.habitat$Habitat, mcm.habitat$Climate) 
chisq.test(tbl)
### because there are many cells with really small values, use below to perform the statistics
ctbl = cbind(tbl[,"Arid"]+tbl[,"Cold"], tbl[,"Temperate"] + tbl[,"Tropical"]) 
chisq.test(ctbl)

# after the file written, the code is now commented
# write.csv(mcm.habitat.proportion,file="Habitat_ProportionData_3April2017.csv")

############################################
##### EDisp and ERich based on global dataset ####
############################################

## Global disparity and diversity
ed.global<-data.frame()
dep.global<-data.frame()
cn<-unique(mcm.habitat$Community.No) ## The number of commuities
ed.global.data<-data.frame()
mean.dep.community<-data.frame()
dep.community.all<-data.frame()
dep.community<-data.frame()
ed.global.data.community<-data.frame()

## disparity parameter of body size
bs.temp.global<-data.frame()
bs.community<-data.frame()
bs.global<-data.frame()
mean.bs.community<-data.frame()

## disparity parameter of diet
dt.temp.global<-data.frame()
dt.community<-data.frame()
dt.global<-data.frame()
mean.dt.community<-data.frame()

## disparity parameter of locomotor mode
lm.temp.global<-data.frame()
lm.community<-data.frame()
lm.global<-data.frame()
mean.lm.community<-data.frame()

for (j in cn) {
  dep.community<-data.frame()
  newdata.global<-mcm.habitat[mcm.habitat$Community.No==j,]
  ed.global.data<-rbind(ed.global.data,c(newdata.global[c(3,4,5)])) ### entire ecological raw dataset
  ed.global.data.community<-rbind(ed.global.data.community,c(newdata.global[c(1,3,4,5)]))
  ### Ecologoical diversity
  data.global<-newdata.global[c(3,4,5)]
  result.global<-nrow(unique(data.global))
  ed.global<-rbind(ed.global,c(result.global))   ## store the data in sequence
  ### Disparity and body size differences
  for (u in 1:(dim(data.global)[1]-1)) {
    for (v in u:(dim(data.global)[1]-1)) {
      ### disparity
      results.temp.global<-abs(data.global[u,]-data.global[v+1,])
      dep.temp.global<-rowSums(results.temp.global)
      dep.community<-rbind(dep.community,c(dep.temp.global))
      dep.global<-rbind(dep.global,c(dep.temp.global))
      ### body size disparity
      bs.temp.global<-abs(data.global[u,1]-data.global[v+1,1])
      bs.community<-rbind(bs.community,c(bs.temp.global))
      bs.global<-rbind(bs.global,c(bs.temp.global))
      ### diet disparity
      dt.temp.global<-abs(data.global[u,2]-data.global[v+1,2])
      dt.community<-rbind(dt.community,c(dt.temp.global))
      dt.global<-rbind(dt.global,c(dt.temp.global))
      ### locomotor mode disparity
      lm.temp.global<-abs(data.global[u,3]-data.global[v+1,3])
      lm.community<-rbind(lm.community,c(lm.temp.global))
      lm.global<-rbind(lm.global,c(lm.temp.global))
    }
  }
  mean.bs.community<-rbind(mean.bs.community,c(mean(bs.community[,1])))
  mean.dt.community<-rbind(mean.dt.community,c(mean(dt.community[,1])))
  mean.lm.community<-rbind(mean.lm.community,c(mean(lm.community[,1])))
  mean.dep.community<-rbind(mean.dep.community, c(mean(dep.community[,1])))
}

mean.dep.global <- mean (dep.global[,1])

### write the EDisp and ERich into same file
EDisp_Erich_together<- cbind(mean.dep.community, ed.global)

## ## Reading the MasterSheet data for the environmental types
Community.proportion.data <- read.csv("EcoHabitat_Proportion_MasterSheet_6June2017.csv", header =T, sep=",")
EDisp_Erich_together_environs <- cbind(EDisp_Erich_together, 
                              Community.proportion.data$Habitat, 
                              Community.proportion.data$Climate, 
                              Community.proportion.data$Vegetation)

colnames(EDisp_Erich_together_environs) <- c("Mean.EDisp", "ERich", "Habitat", "Climate", "Vegetation")

## ggplot EDisp vs ERich with annatations of different environmental types

## no annotation
EDisp_ERich_plot <- ggplot(EDisp_Erich_together_environs, aes(Mean.EDisp, ERich)) +
  geom_point(color="grey20",fill="black", size = 3)

## two habitats
EDisp_ERich_habitat <- ggplot(EDisp_Erich_together_environs, aes(Mean.EDisp, ERich)) + 
  geom_point(aes(color=Habitat, shape=Habitat), size = 3) + 
  scale_color_manual(values = c("grey30","black")) +  ## grey30 is closed habitat
  scale_shape_manual(values = c(4, 15))               ## 4 is closed habitat


## four climates
EDisp_ERich_climate <- ggplot(EDisp_Erich_together_environs, aes(Mean.EDisp, ERich)) + 
  geom_point(aes(color=Climate, shape=Climate), size = 3) + 
  scale_color_manual(values = c("grey20",         ### Arid; 
                                "dodgerblue3",    ## cold;
                                "green3",         ## temperate;
                                "darkorange")) +  ## tropical
  scale_shape_manual(values = c(1, 15, 4, 5))    ## same sequence


## eight vegetations
EDisp_ERich_vegetation <- ggplot(EDisp_Erich_together_environs, aes(Mean.EDisp, ERich)) + 
  geom_point(aes(color=Vegetation, shape=Vegetation), size = 3) + 
  scale_color_manual(values = c("darkviolet",       # boreal forest
                                "khaki",       # desert
                                "goldenrod",     # grassland
                                "firebrick1",         # savanna
                                "darkgoldenrod1",       # shrubland
                                "green3",        # temperate forest 
                                "red",         # tropical rainforest
                                "orangered")) + # tropical seasonal forest
  scale_shape_manual(values = c(1,        # boreal forest
                                18,        # desert
                                15,     # grassland
                                2,         # savanna
                                4,      # shrubland
                                0,       # temperate forest 
                                16,        # tropical rainforest
                                11))     # tropical seasonal forest
grid.arrange(EDisp_ERich_plot , EDisp_ERich_climate,
             EDisp_ERich_habitat, EDisp_ERich_vegetation, 
             ncol=2, nrow=2)

######################################################################
## relationships between EDisp and ERich in different environments
#####################################################################

## data of different environments
EDisp_ERich_open <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Habitat=="Open",]
EDisp_ERich_close <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Habitat=="Close",]

EDisp_ERich_tropical <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Climate=="Tropical",]
EDisp_ERich_arid <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Climate=="Arid",]
EDisp_ERich_temperate <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Climate=="Temperate",]
EDisp_ERich_cold <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Climate=="Cold",]

EDisp_ERich_tropical_rainforest <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Tropical rainforest",]
EDisp_ERich_tropical_seasonal_forest <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Tropical seasonal forest",]
EDisp_ERich_savanna <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Savanna",]
EDisp_ERich_grassland <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Grassland",]
EDisp_ERich_shrubland <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Shrubland",]
EDisp_ERich_desert <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Desert",]
EDisp_ERich_temperate_forest <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Temperate forest",]
EDisp_ERich_boreal_forest <- EDisp_Erich_together_environs[EDisp_Erich_together_environs$Vegetation=="Boreal forest",]

## Regressions of above
Regression_results_open <- summary(lm(EDisp_ERich_open$ERich~EDisp_ERich_open$Mean.EDisp))
Regression_results_close <- summary(lm(EDisp_ERich_close$ERich~EDisp_ERich_close$Mean.EDisp))

Regression_results_tropical <- summary(lm(EDisp_ERich_tropical$ERich~EDisp_ERich_tropical$Mean.EDisp))
Regression_results_arid <- summary(lm(EDisp_ERich_arid$ERich~EDisp_ERich_arid$Mean.EDisp))
Regression_results_temperate <- summary(lm(EDisp_ERich_temperate$ERich~EDisp_ERich_temperate$Mean.EDisp))
Regression_results_cold <- summary(lm(EDisp_ERich_cold$ERich~EDisp_ERich_cold$Mean.EDisp))

Regression_results_tropical_rainforest <- summary(lm(EDisp_ERich_tropical_rainforest$ERich~EDisp_ERich_tropical_rainforest$Mean.EDisp))
Regression_results_tropical_seasonal_forest <- summary(lm(EDisp_ERich_tropical_seasonal_forest$ERich~EDisp_ERich_tropical_seasonal_forest$Mean.EDisp))
Regression_results_savanna <- summary(lm(EDisp_ERich_savanna$ERich~EDisp_ERich_savanna$Mean.EDisp))
Regression_results_grassland <- summary(lm(EDisp_ERich_grassland$ERich~EDisp_ERich_grassland$Mean.EDisp))
Regression_results_shrubland<- summary(lm(EDisp_ERich_shrubland$ERich~EDisp_ERich_shrubland$Mean.EDisp))
Regression_results_desert <- summary(lm(EDisp_ERich_desert$ERich~EDisp_ERich_desert$Mean.EDisp))
Regression_results_temperate_forest <- summary(lm(EDisp_ERich_temperate_forest$ERich~EDisp_ERich_temperate_forest$Mean.EDisp))
Regression_results_boreal_forest <- summary(lm(EDisp_ERich_boreal_forest$ERich~EDisp_ERich_boreal_forest$Mean.EDisp))

## function of the stat of regression
Regression_stat_func <- function(input) {
  stat.table <- cbind(input$r.squared, input$fstatistic[1], input$coefficients[2,4])
  return(stat.table)
}

Regression_results_stat_all <- rbind(Regression_stat_func(Regression_results_close),
                                     Regression_stat_func(Regression_results_open),
                                     Regression_stat_func(Regression_results_tropical),
                                     Regression_stat_func(Regression_results_arid),
                                     Regression_stat_func(Regression_results_temperate),
                                     Regression_stat_func(Regression_results_cold),
                                     Regression_stat_func(Regression_results_tropical_rainforest),
                                     Regression_stat_func(Regression_results_tropical_seasonal_forest),
                                     Regression_stat_func(Regression_results_savanna),
                                     Regression_stat_func(Regression_results_grassland),
                                     Regression_stat_func(Regression_results_shrubland),
                                     Regression_stat_func(Regression_results_desert),
                                     Regression_stat_func(Regression_results_temperate_forest),
                                     Regression_stat_func(Regression_results_boreal_forest))

colnames(Regression_results_stat_all) <- c("R-squared", "F", "P")
rownames(Regression_results_stat_all) <- c("Close", "Open",
                                           "Tropical", "Arid", "Temperate", "Cold",
                                           "Tropical rain forest", "Tropical seasonal forest", 
                                           "Savanna", "Grassland", "Shrubland", "Desert",
                                           "Temperate forest", "Boreal forest")

write.csv(Regression_results_stat_all, file = "EDisp_ERich_Regression_stat_table_16Oct2017.csv")

#write.csv(EDisp_Erich_together, file = "EDisp_ERich_30June2017.csv")


################################################################################################
##### kendall's concordance analysis for independence
################################################################################################
kendall.global(ed.global.data[,1:3], nperm = 999, mult = "holm")
kendall.post(ed.global.data[,1:3], nperm = 999, mult = "holm")
cor(ed.global.data[,1:3], method="kendall", use="pairwise")


#########################################################
##### EDisp + ERich Aanalyses based on two different habitat types ####
#########################################################

#### open habitat
mcm.open<-mcm.habitat[mcm.habitat$Habitat == "Open",]
write.table(mcm.open[,3:5], file="data.open.csv", sep=",", row.names=F, col.names = F)
#### close habitat
mcm.close<-mcm.habitat[mcm.habitat$Habitat == "Close",]
write.table(mcm.close[,3:5], file="data.close.csv", sep=",", row.names=F, col.names = F)

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
      ## body size disparity
      bs.temp.open<-abs(data.open[u.open,1]-data.open[v.open+1,1])
      bs.community.open<-rbind(bs.community.open,c(bs.temp.open))
      bs.open<-rbind(bs.open,c(bs.temp.open))
      ### diet disparity
      dt.temp.open<-abs(data.open[u.open,2]-data.open[v.open+1,2])
      dt.community.open<-rbind(dt.community.open,c(dt.temp.open))
      dt.open<-rbind(dt.open,c(dt.temp.open))
      ### locomotor mode disparity
      lm.temp.open<-abs(data.open[u.open,3]-data.open[v.open+1,3])
      lm.community.open<-rbind(lm.community.open,c(lm.temp.open))
      lm.open<-rbind(lm.open,c(lm.temp.open))
    }
  }
  mean.bs.community.open<-rbind(mean.bs.community.open,c(mean(bs.community.open[,1])))
  mean.dt.community.open<-rbind(mean.dt.community.open,c(mean(dt.community.open[,1])))
  mean.lm.community.open<-rbind(mean.lm.community.open,c(mean(lm.community.open[,1])))
}

ed.open.data$Eco.type <- with(ed.open.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

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
      ## body size disparity
      bs.temp.close<-abs(data.close[u.close,1]-data.close[v.close+1,1])
      bs.community.close<-rbind(bs.community.close,c(bs.temp.close))
      bs.close<-rbind(bs.close,c(bs.temp.close))
      ### diet disparity
      dt.temp.close<-abs(data.close[u.close,2]-data.close[v.close+1,2])
      dt.community.close<-rbind(dt.community.close,c(dt.temp.close))
      dt.close<-rbind(dt.close,c(dt.temp.close))
      ### locomotor mode disparity
      lm.temp.close<-abs(data.close[u.close,3]-data.close[v.close+1,3])
      lm.community.close<-rbind(lm.community.close,c(lm.temp.close))
      lm.close<-rbind(lm.close,c(lm.temp.close))
    }
  }
  mean.bs.community.close<-rbind(mean.bs.community.close,c(mean(bs.community.close[,1])))
  mean.dt.community.close<-rbind(mean.dt.community.close,c(mean(dt.community.close[,1])))
  mean.lm.community.close<-rbind(mean.lm.community.close,c(mean(lm.community.close[,1])))
}

ed.close.data$Eco.type <- with(ed.close.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

combined.dep.close<-matrix(c(bs.close[,1], dt.close[,1], lm.close[,1]), ncol=3, byrow = F)
combined.dep.open<-matrix(c(bs.open[,1], dt.open[,1], lm.open[,1]), ncol=3, byrow = F)

### Separate plots

par(mfcol=c(3,2))

for (tr in 1:dim(combined.dep.close)[2]) {
  barplot(table(combined.dep.close[,tr])/dim(combined.dep.close)[1], xlim=c(0,9))
}

for (tr in 1:dim(combined.dep.open)[2]) {
  barplot(table(combined.dep.open[,tr])/dim(combined.dep.open)[1], xlim=c(0,9))
} 

### Combined plots

par(mfcol=c(2,1))

bs.disp.close<-table(combined.dep.close[,1])
dt.disp.close<-table(combined.dep.close[,2])
lm.disp.close<-table(combined.dep.close[,3])
barplot(rbind(bs.disp.close,dt.disp.close,lm.disp.close)/dim(combined.dep.close)[1], col=c("red","orange","blue"), space=1)


bs.disp.open<-table(combined.dep.open[,1])
dt.disp.open<-table(combined.dep.open[,2])
lm.disp.open<-table(combined.dep.open[,3])
barplot(rbind(bs.disp.open,dt.disp.open,lm.disp.open)/dim(combined.dep.open)[1], col=c("red","orange","blue"), space=1)


#### disparity contribution of each eoclogical paratmer
par(mfcol=c(2,1))
barplot(c(sum(combined.dep.close[,1]),sum(combined.dep.close[,2]),sum(combined.dep.close[,3]))/sum(combined.dep.close))
barplot(c(sum(combined.dep.open[,1]),sum(combined.dep.open[,2]),sum(combined.dep.open[,3]))/sum(combined.dep.open))

#### stack plots
close.type<-rbind(sum(combined.dep.close[,1])/sum(combined.dep.close), sum(combined.dep.close[,2])/sum(combined.dep.close), 
          sum(combined.dep.close[,3])/sum(combined.dep.close))
open.type<-rbind(sum(combined.dep.open[,1])/sum(combined.dep.open), sum(combined.dep.open[,2])/sum(combined.dep.open), 
          sum(combined.dep.open[,3])/sum(combined.dep.open))
combined.type<-cbind(close.type*mean(dep.close[,1]),open.type*mean(dep.open[,1]))
barplot(combined.type, space=1, col=c("dodgerblue3","green3","darkorange"))

### Proportional dataset

par(mfcol=c(2,3),
    oma=c(1,1,0,1), 
    mar=c(2,2,1,0), 
    mgp=c(1,0.6,0),
    tck=-0.02)


# body size
barplot(table(ed.close.data[,1])/sum(table(ed.close.data[,1])), space=1, ylim=c(0,0.6), col="grey60")              # close
barplot(table(ed.open.data[,1])/sum(table(ed.open.data[,1])), space=1, ylim=c(0,0.6), col="grey30")                # open

# diet
barplot(table(ed.close.data[,2])/sum(table(ed.close.data[,2])), space=1, ylim=c(0,0.6), col="grey60")              # close
barplot(table(ed.open.data[,2])/sum(table(ed.open.data[,2])), space=1, ylim=c(0,0.6), col="grey30")                # open

# locomotor mode
barplot(table(ed.close.data[,3])/sum(table(ed.close.data[,3])), space=1, ylim=c(0,0.6), col="grey60")              # close
barplot(table(ed.open.data[,3])/sum(table(ed.open.data[,3])), space=1, ylim=c(0,0.6), col="grey30")                # open


### close vs open on the EDisp plots
par(mfrow=c(1,1))
error.bar <- function(x, y, upper, lower=upper, length=0.1,...) {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

mean.habitat <- cbind(mean(dep.close[,1]), mean(dep.open[,1]))
sd.habitat <- cbind(sd(dep.close[,1]),sd(dep.open[,1]))
mean.habitat.names <- c("Close","Open")
length.habitat <- cbind(length(t(dep.close)),length(t(dep.open)))

barx <- barplot(mean.habitat, names.arg=mean.habitat.names,ylim=c(0,6), col="gold3", 
                cex.axis=0.8,axis.lty=1, space=1, cex.names=0.8,las=1)
## plot the error bar
error.bar(barx, mean.habitat, 1.96*sd.habitat/sqrt(length.habitat))

##### dissimilarity between close and open habitats
dis.close.open<-1-length(intersect(unique(ed.close.data$Eco.type), unique(ed.open.data$Eco.type)))/length(union(unique(ed.close.data$Eco.type), unique(ed.open.data$Eco.type)))

dis.close.open

###############################################################
##### EDisp + Erich based on four different climatic types ####
###############################################################

mcm.tropical.climate<-mcm.habitat[mcm.habitat$Climate == "Tropical",]
write.table(mcm.tropical.climate[,3:5], file="data.tropical.climate.csv", sep=",", row.names = F, col.names = F)

mcm.arid.climate<-mcm.habitat[mcm.habitat$Climate == "Arid",]
write.table(mcm.arid.climate[,3:5], file="data.arid.climate.csv", sep=",", row.names = F, col.names = F)

mcm.temperate.climate<-mcm.habitat[mcm.habitat$Climate== "Temperate",]
write.table(mcm.temperate.climate[,3:5], file="data.temperate.climate.csv", sep=",", row.names = F, col.names = F)

mcm.cold.climate<-mcm.habitat[mcm.habitat$Climate == "Cold",]
write.table(mcm.cold.climate[,3:5], file="data.cold.climate.csv", sep=",", row.names = F, col.names = F)

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
      ## body size disparity
      bs.temp.tropical.climate<-abs(data.tropical.climate[u.tropical.climate,1]-data.tropical.climate[v.tropical.climate+1,1])
      bs.community.tropical.climate<-rbind(bs.community.tropical.climate,c(bs.temp.tropical.climate))
      bs.tropical.climate<-rbind(bs.tropical.climate,c(bs.temp.tropical.climate))
      ### diet disparity
      dt.temp.tropical.climate<-abs(data.tropical.climate[u.tropical.climate,2]-data.tropical.climate[v.tropical.climate+1,2])
      dt.community.tropical.climate<-rbind(dt.community.tropical.climate,c(dt.temp.tropical.climate))
      dt.tropical.climate<-rbind(dt.tropical.climate,c(dt.temp.tropical.climate))
      ### locomotor mode disparity
      lm.temp.tropical.climate<-abs(data.tropical.climate[u.tropical.climate,3]-data.tropical.climate[v.tropical.climate+1,3])
      lm.community.tropical.climate<-rbind(lm.community.tropical.climate,c(lm.temp.tropical.climate))
      lm.tropical.climate<-rbind(lm.tropical.climate,c(lm.temp.tropical.climate))
    }
  }
  mean.bs.community.tropical.climate<-rbind(mean.bs.community.tropical.climate,c(mean(bs.community.tropical.climate[,1])))
  mean.dt.community.tropical.climate<-rbind(mean.dt.community.tropical.climate,c(mean(dt.community.tropical.climate[,1])))
  mean.lm.community.tropical.climate<-rbind(mean.lm.community.tropical.climate,c(mean(lm.community.tropical.climate[,1])))
}

ed.tropical.climate.data$Eco.type <- with(ed.tropical.climate.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))


## combining bs, dt, and lm
combined.dep.tropical.climate<-matrix(c(bs.tropical.climate[,1], dt.tropical.climate[,1], lm.tropical.climate[,1]), ncol=3, byrow = F)

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
      ## body size disparity
      bs.temp.arid.climate<-abs(data.arid.climate[u.arid.climate,1]-data.arid.climate[v.arid.climate+1,1])
      bs.community.arid.climate<-rbind(bs.community.arid.climate,c(bs.temp.arid.climate))
      bs.arid.climate<-rbind(bs.arid.climate,c(bs.temp.arid.climate))
      ### diet disparity
      dt.temp.arid.climate<-abs(data.arid.climate[u.arid.climate,2]-data.arid.climate[v.arid.climate+1,2])
      dt.community.arid.climate<-rbind(dt.community.arid.climate,c(dt.temp.arid.climate))
      dt.arid.climate<-rbind(dt.arid.climate,c(dt.temp.arid.climate))
      ### locomotor mode disparity
      lm.temp.arid.climate<-abs(data.arid.climate[u.arid.climate,3]-data.arid.climate[v.arid.climate+1,3])
      lm.community.arid.climate<-rbind(lm.community.arid.climate,c(lm.temp.arid.climate))
      lm.arid.climate<-rbind(lm.arid.climate,c(lm.temp.arid.climate))
    }
  }
  mean.bs.community.arid.climate<-rbind(mean.bs.community.arid.climate,c(mean(bs.community.arid.climate[,1])))
  mean.dt.community.arid.climate<-rbind(mean.dt.community.arid.climate,c(mean(dt.community.arid.climate[,1])))
  mean.lm.community.arid.climate<-rbind(mean.lm.community.arid.climate,c(mean(lm.community.arid.climate[,1])))
}

ed.arid.climate.data$Eco.type <- with(ed.arid.climate.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

## combining bs, dt, and lm
combined.dep.arid.climate<-matrix(c(bs.arid.climate[,1], dt.arid.climate[,1], lm.arid.climate[,1]), ncol=3, byrow = F)

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
      ## body size disparity
      bs.temp.temperate.climate<-abs(data.temperate.climate[u.temperate.climate,1]-data.temperate.climate[v.temperate.climate+1,1])
      bs.community.temperate.climate<-rbind(bs.community.temperate.climate,c(bs.temp.temperate.climate))
      bs.temperate.climate<-rbind(bs.temperate.climate,c(bs.temp.temperate.climate))
      ### diet disparity
      dt.temp.temperate.climate<-abs(data.temperate.climate[u.temperate.climate,2]-data.temperate.climate[v.temperate.climate+1,2])
      dt.community.temperate.climate<-rbind(dt.community.temperate.climate,c(dt.temp.temperate.climate))
      dt.temperate.climate<-rbind(dt.temperate.climate,c(dt.temp.temperate.climate))
      ### locomotor mode disparity
      lm.temp.temperate.climate<-abs(data.temperate.climate[u.temperate.climate,3]-data.temperate.climate[v.temperate.climate+1,3])
      lm.community.temperate.climate<-rbind(lm.community.temperate.climate,c(lm.temp.temperate.climate))
      lm.temperate.climate<-rbind(lm.temperate.climate,c(lm.temp.temperate.climate))
    }
  }
  mean.bs.community.temperate.climate<-rbind(mean.bs.community.temperate.climate,c(mean(bs.community.temperate.climate[,1])))
  mean.dt.community.temperate.climate<-rbind(mean.dt.community.temperate.climate,c(mean(dt.community.temperate.climate[,1])))
  mean.lm.community.temperate.climate<-rbind(mean.lm.community.temperate.climate,c(mean(lm.community.temperate.climate[,1])))
}

ed.temperate.climate.data$Eco.type <- with(ed.temperate.climate.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

## combining bs, dt, and lm
combined.dep.temperate.climate<-matrix(c(bs.temperate.climate[,1], dt.temperate.climate[,1], lm.temperate.climate[,1]), ncol=3, byrow = F)

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
      ## body size disparity
      bs.temp.cold.climate<-abs(data.cold.climate[u.cold.climate,1]-data.cold.climate[v.cold.climate+1,1])
      bs.community.cold.climate<-rbind(bs.community.cold.climate,c(bs.temp.cold.climate))
      bs.cold.climate<-rbind(bs.cold.climate,c(bs.temp.cold.climate))
      ### diet disparity
      dt.temp.cold.climate<-abs(data.cold.climate[u.cold.climate,2]-data.cold.climate[v.cold.climate+1,2])
      dt.community.cold.climate<-rbind(dt.community.cold.climate,c(dt.temp.cold.climate))
      dt.cold.climate<-rbind(dt.cold.climate,c(dt.temp.cold.climate))
      ### locomotor mode disparity
      lm.temp.cold.climate<-abs(data.cold.climate[u.cold.climate,3]-data.cold.climate[v.cold.climate+1,3])
      lm.community.cold.climate<-rbind(lm.community.cold.climate,c(lm.temp.cold.climate))
      lm.cold.climate<-rbind(lm.cold.climate,c(lm.temp.cold.climate))
    }
  }
  mean.bs.community.cold.climate<-rbind(mean.bs.community.cold.climate,c(mean(bs.community.cold.climate[,1])))
  mean.dt.community.cold.climate<-rbind(mean.dt.community.cold.climate,c(mean(dt.community.cold.climate[,1])))
  mean.lm.community.cold.climate<-rbind(mean.lm.community.cold.climate,c(mean(lm.community.cold.climate[,1])))
}

ed.cold.climate.data$Eco.type <- with(ed.cold.climate.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

## combining bs, dt, and lm
combined.dep.cold.climate<-matrix(c(bs.cold.climate[,1], dt.cold.climate[,1], lm.cold.climate[,1]), ncol=3, byrow = F)

### save the mean and sd of EDisp and ERich
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

write.csv(climate.dep, file = "EDisp_climate_statistic_24July2017.csv")
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

write.csv(climate.ed, file = "ERich_climate_statistic_24July2017.csv")


## See the summary
eco.type.tropical.climate <- table(ed.tropical.climate.data$Eco.type)
eco.type.arid.climate <- table(ed.arid.climate.data$Eco.type)
eco.type.temperate.climate <- table(ed.temperate.climate.data$Eco.type)
eco.type.cold.climate <- table(ed.cold.climate.data$Eco.type)


#### disparity contribution of each eoclogical paratmer
par(mfcol=c(2,2))

#### stack plots
tropical.climatic.type<-rbind(sum(bs.tropical.climate[,1])/sum(combined.dep.tropical.climate), sum(dt.tropical.climate[,1])/sum(combined.dep.tropical.climate), 
                          sum(lm.tropical.climate[,1])/sum(combined.dep.tropical.climate))
arid.climatic.type<-rbind(sum(bs.arid.climate[,1])/sum(combined.dep.arid.climate), sum(dt.arid.climate[,1])/sum(combined.dep.arid.climate), 
                          sum(lm.arid.climate[,1])/sum(combined.dep.arid.climate))
temperate.climatic.type<-rbind(sum(bs.temperate.climate[,1])/sum(combined.dep.temperate.climate), sum(dt.temperate.climate[,1])/sum(combined.dep.temperate.climate), 
                          sum(lm.temperate.climate[,1])/sum(combined.dep.temperate.climate))
cold.climatic.type<-rbind(sum(bs.cold.climate[,1])/sum(combined.dep.cold.climate), sum(dt.cold.climate[,1])/sum(combined.dep.cold.climate), 
                  sum(lm.cold.climate[,1])/sum(combined.dep.cold.climate))

#### stack plots
combined.climatic.type<-cbind(tropical.climatic.type*mean(dep.tropical.climate[,1]),
                    arid.climatic.type*mean(dep.arid.climate[,1]),
                    temperate.climatic.type*mean(dep.temperate.climate[,1]),
                    cold.climatic.type*mean(dep.cold.climate[,1]))
barplot(combined.climatic.type, space=1, col=c("dodgerblue3","green3","darkorange"))


# ################################################
# ## plot of the each value of the functional trait
# #################################################
# par(mfcol=c(4,3),
#     oma=c(1,1,0,1), 
#     mar=c(2,4,1,1), 
#     mgp=c(1,0.6,0),
#     tck=-0.02)
# 
# # body size
# barplot(table(ed.tropical.climate.data[,1]), ylim=c(0,130), space=1, col="darkorange") # tropical
# barplot(table(ed.arid.climate.data[,1]), space=1, ylim=c(0,130), col="grey20")              # arid
# barplot(table(ed.temperate.climate.data[,1]), space=1, ylim=c(0,130), col="green3")         # temperate
# barplot(table(ed.cold.climate.data[,1]), space=1, ylim=c(0,130), col="dodgerblue3")              # cold
# 
# # diet
# barplot(table(ed.tropical.climate.data[,2]),space=1,ylim=c(0,130), col="darkorange") # tropical
# barplot(table(ed.arid.climate.data[,2]),space=1, ylim=c(0,130), col="grey20")              # arid
# barplot(table(ed.temperate.climate.data[,2]),space=1, ylim=c(0,130), col="green3")         # temperate
# barplot(table(ed.cold.climate.data[,2]),space=1, ylim=c(0,130), col="dodgerblue3")              # cold
# 
# # locomotor mode
# barplot(table(ed.tropical.climate.data[,3]),space=1, ylim=c(0,130), col="darkorange") # tropical
# barplot(table(ed.arid.climate.data[,3]),space=1, ylim=c(0,130), col="grey20")              # arid
# barplot(table(ed.temperate.climate.data[,3]),space=1, ylim=c(0,130), col="green3")         # temperate
# barplot(table(ed.cold.climate.data[,3]),space=1, ylim=c(0,130), col="dodgerblue3")              # cold

### Plots of proportoin data of four climatic types

par(mfcol=c(4,3),
    oma=c(1,1,0,1), 
    mar=c(2,2,1,0), 
    mgp=c(1,0.6,0),
    tck=-0.02)

# body size
barplot(table(ed.tropical.climate.data[,1])/sum(table(ed.tropical.climate.data[,1])), ylim=c(0,0.6), space=1, col="darkorange") # tropical
barplot(table(ed.arid.climate.data[,1])/sum(table(ed.arid.climate.data[,1])), space=1, ylim=c(0,0.6), col="grey50")              # arid
barplot(table(ed.temperate.climate.data[,1])/sum(table(ed.temperate.climate.data[,1])), space=1, ylim=c(0,0.6), col="green3")         # temperate
barplot(table(ed.cold.climate.data[,1])/sum(table(ed.cold.climate.data[,1])), space=1, ylim=c(0,0.6), col="dodgerblue3")              # cold

# diet
barplot(table(ed.tropical.climate.data[,2])/sum(table(ed.tropical.climate.data[,2])),space=1,ylim=c(0,0.6), col="darkorange") # tropical
barplot(table(ed.arid.climate.data[,2])/sum(table(ed.arid.climate.data[,2])),space=1, ylim=c(0,0.6), col="grey50")              # arid
barplot(table(ed.temperate.climate.data[,2])/sum(table(ed.temperate.climate.data[,2])),space=1, ylim=c(0,0.6), col="green3")         # temperate
barplot(table(ed.cold.climate.data[,2])/sum(table(ed.cold.climate.data[,2])),space=1, ylim=c(0,0.6), col="dodgerblue3")              # cold

# locomotor mode
barplot(table(ed.tropical.climate.data[,3])/sum(table(ed.tropical.climate.data[,3])),space=1, ylim=c(0,0.6), col="darkorange") # tropical
barplot(table(ed.arid.climate.data[,3])/sum(table(ed.arid.climate.data[,3])),space=1, ylim=c(0,0.6), col="grey50")              # arid
barplot(table(ed.temperate.climate.data[,3])/sum(table(ed.temperate.climate.data[,3])),space=1, ylim=c(0,0.6), col="green3")         # temperate
barplot(table(ed.cold.climate.data[,3])/sum(table(ed.cold.climate.data[,3])),space=1, ylim=c(0,0.6), col="dodgerblue3")              # cold

##### dissimilarity among each climate
dis.tropical.arid<-1-length(intersect(unique(ed.tropical.climate.data$Eco.type), unique(ed.arid.climate.data$Eco.type)))/length(union(unique(ed.tropical.climate.data$Eco.type), unique(ed.arid.climate.data$Eco.type)))

dis.tropical.temperate<-1-length(intersect(unique(ed.tropical.climate.data$Eco.type), unique(ed.temperate.climate.data$Eco.type)))/length(union(unique(ed.tropical.climate.data$Eco.type), unique(ed.temperate.climate.data$Eco.type)))

dis.tropical.cold<-1-length(intersect(unique(ed.tropical.climate.data$Eco.type), unique(ed.cold.climate.data$Eco.type)))/length(union(unique(ed.tropical.climate.data$Eco.type), unique(ed.cold.climate.data$Eco.type)))

dis.arid.temperate<-1-length(intersect(unique(ed.temperate.climate.data$Eco.type), unique(ed.arid.climate.data$Eco.type)))/length(union(unique(ed.temperate.climate.data$Eco.type),unique(ed.arid.climate.data$Eco.type)))

dis.arid.cold<-1-length(intersect(unique(ed.cold.climate.data$Eco.type), unique(ed.arid.climate.data$Eco.type)))/length(union(unique(ed.cold.climate.data$Eco.type), unique(ed.arid.climate.data$Eco.type)))

dis.temperate.cold<-1-length(intersect(unique(ed.cold.climate.data$Eco.type), unique(ed.temperate.climate.data$Eco.type)))/length(union(unique(ed.cold.climate.data$Eco.type),unique(ed.temperate.climate.data$Eco.type)))

dis.tropical.arid
dis.tropical.temperate
dis.tropical.cold
dis.arid.temperate
dis.arid.cold
dis.temperate.cold


################################
# MEAN DISPARITY by different climates

## tropical.climate Mean Disparity
mean.tropical.climate<-data.frame()
sd.tropical.climate<-data.frame()

mean.temp.tropical.climate<-apply(as.matrix(dep.tropical.climate[,1]),2,mean)
mean.tropical.climate<-rbind(mean.tropical.climate,c(mean.temp.tropical.climate))
sd.temp.tropical.climate<-apply(as.matrix(dep.tropical.climate[,1]),2,sd)
sd.tropical.climate<-rbind(sd.tropical.climate,c(sd.temp.tropical.climate))

mean.tropical.climate<-t(mean.tropical.climate)
sd.tropical.climate<-t(sd.tropical.climate)

## arid.climate Mean Disparity
mean.arid.climate<-data.frame()
sd.arid.climate<-data.frame()

mean.temp.arid.climate<-apply(as.matrix(dep.arid.climate[,1]),2,mean)
mean.arid.climate<-rbind(mean.arid.climate,c(mean.temp.arid.climate))
sd.temp.arid.climate<-apply(as.matrix(dep.arid.climate[,1]),2,sd)
sd.arid.climate<-rbind(sd.arid.climate,c(sd.temp.arid.climate))

mean.arid.climate<-t(mean.arid.climate)
sd.arid.climate<-t(sd.arid.climate)

## temperate.climate Mean Disparity
mean.temperate.climate<-data.frame()
sd.temperate.climate<-data.frame()

mean.temp.temperate.climate<-apply(as.matrix(dep.temperate.climate[,1]),2,mean)
mean.temperate.climate<-rbind(mean.temperate.climate,c(mean.temp.temperate.climate))
sd.temp.temperate.climate<-apply(as.matrix(dep.temperate.climate[,1]),2,sd)
sd.temperate.climate<-rbind(sd.temperate.climate,c(sd.temp.temperate.climate))

mean.temperate.climate<-t(mean.temperate.climate)
sd.temperate.climate<-t(sd.temperate.climate)

## Cold Mean Disparity
mean.cold.climate<-data.frame()
sd.cold.climate<-data.frame()

mean.temp.cold.climate<-apply(as.matrix(dep.cold.climate[,1]),2,mean)
mean.cold.climate<-rbind(mean.cold.climate,c(mean.temp.cold.climate))
sd.temp.cold.climate<-apply(as.matrix(dep.cold.climate[,1]),2,sd)
sd.cold.climate<-rbind(sd.cold.climate,c(sd.temp.cold.climate))

mean.cold.climate<-t(mean.cold.climate)
sd.cold.climate<-t(sd.cold.climate)

## Barplot of Disparity
## funciton of standard error in barplot
error.bar <- function(x, y, upper, lower=upper, length=0.1,...) {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

mean.all <- cbind(mean.tropical.climate,mean.arid.climate,mean.temperate.climate,mean.cold.climate)
sd.all <- cbind(sd.tropical.climate,sd.arid.climate,sd.temperate.climate,sd.cold.climate)
mean.names <- c("Tropical","Arid","Temperate","Cold")
length.all <- cbind(length(t(dep.tropical.climate)),length(t(dep.arid.climate)),length(t(dep.temperate.climate)),
                    length(t(dep.cold.climate)))

barx <- barplot(mean.all, names.arg=mean.names,ylim=c(0,6), col="grey70", 
                cex.axis=0.8,axis.lty=1, space=1, cex.names=0.8,las=1)
## plot the error bar
error.bar(barx, mean.all, 1.96*sd.all/sqrt(length.all))

mtext(outer=F, side=2, line=1.5, text="Mean Disparity",cex=0.8)
mtext(outer=F, side=3, line=-2, text="Ecological disparity", cex=1.5)


### Ecological diversity
## Tropical Mean Diversity
mean.tropical.climate.ed<-mean(ed.tropical.climate[,1])
sd.tropical.climate.ed<-sd(ed.tropical.climate[,1])
mean.tropical.climate.ed<-t(mean.tropical.climate.ed)
sd.tropical.climate.ed<-t(sd.tropical.climate.ed)

## arid.climate Mean Diversity

mean.arid.climate.ed<-mean(ed.arid.climate[,1])
sd.arid.climate.ed<-sd(ed.arid.climate[,1])
mean.arid.climate.ed<-t(mean.arid.climate.ed)
sd.arid.climate.ed<-t(sd.arid.climate.ed)

## temperate.climate Mean Diversity

mean.temperate.climate.ed<-mean(ed.temperate.climate[,1])
sd.temperate.climate.ed<-sd(ed.temperate.climate[,1])
mean.temperate.climate.ed<-t(mean.temperate.climate.ed)
sd.temperate.climate.ed<-t(sd.temperate.climate.ed)

## cold.climate Mean Diversity

mean.cold.climate.ed<-mean(ed.cold.climate[,1])
sd.cold.climate.ed<-sd(ed.cold.climate[,1])
mean.cold.climate.ed<-t(mean.cold.climate.ed)
sd.cold.climate.ed<-t(sd.cold.climate.ed)


###############################
### plot Diversity together ###
###############################
mean.ed.all<-cbind(mean.tropical.climate.ed[,1],
                   mean.arid.climate.ed[,1],
                   mean.temperate.climate.ed[,1],
                   mean.cold.climate.ed[,1])

sd.ed.all<-cbind(sd.tropical.climate.ed[,1],
                 sd.arid.climate.ed[,1],
                 sd.temperate.climate.ed[,1],
                 sd.cold.climate.ed[,1])

length.ed.all<-t(matrix(cbind(length(ed.tropical.climate[,1]),
                              length(ed.arid.climate[,1]),length(ed.temperate.climate[,1]),
                              length(ed.cold.climate[,1]))))

## all together
barx <- barplot(mean.ed.all, names.arg=c("tropical.climate",
                                         "arid.climate","temperate.climate","cold.climate"), 
                ylim=c(0,15), col="grey70", cex.axis=0.8,
                axis.lty=1, cex.names=0.8,
                las=1,beside=TRUE)
## plot the error bar
error.bar(barx, mean.ed.all, 1.96*sd.ed.all/sqrt(length.ed.all))
mtext(outer=F, side=2, line=1.5, text="Mean Diversity",cex=0.8)


##################################################
# Barplots of ecological disparity and diversity 
##################################################

par(mfcol=c(1,1),
    oma=c(1,1,0,1), 
    mar=c(2,4,0,1), 
    mgp=c(1,0.6,0),
    tck=-0.02)

mean.together<-rbind(mean.all, mean.ed.all)
sd.together<-rbind(sd.all,sd.ed.all)
length.together<-rbind(length.all, length.ed.all)
barx <- barplot(mean.together, beside=T,names.arg=c("Tropical",
                                                    "Arid","Temperate","Cold"),
                ylim=c(0,15), col=c("gold3","cadetblue"), cex.axis=1.2,
                axis.lty=1, cex.names=1.5,las=1, width = 0.75, xlim = c(0,10))
## plot the error bar
error.bar(barx, mean.together, 1.96*sd.together/sqrt(length.together))
mtext(outer=F, side=2, line=2, text="Mean Value",cex=1.5)
legend(0, 16,legend=c("Ecological disparity", "Ecological diversity"),
       fill=c("gold3","cadetblue"), horiz=T,cex=1.2,box.col = "white", border="white")


####################################################
##### EDisp + Erich based on eight vegetation types ####
####################################################

##### The data of each habitat
# Boreal forest
mcm.boreal.forest<-mcm.habitat[mcm.habitat$Vegetation == "Boreal forest",]
write.table(mcm.boreal.forest[,3:5], file="data.boreal.forest.csv", sep=",", row.names=F, col.names = F)
# Desert
mcm.desert <-mcm.habitat[mcm.habitat$Vegetation == "Desert",]
write.table(mcm.desert[,3:5], file="data.desert.csv", sep=",", row.names=F, col.names = F)
# Temperate forest
mcm.temperate.forest<-mcm.habitat[mcm.habitat$Vegetation == "Temperate forest",]
write.table(mcm.temperate.forest[,3:5], file="data.temperate.forest.csv", sep=",", row.names=F, col.names = F)
# Grassland
mcm.grassland<-mcm.habitat[mcm.habitat$Vegetation== "Grassland",]
write.table(mcm.grassland[,3:5], file="data.grassland.csv", sep=",", row.names=F, col.names = F)
# Shrubland
mcm.shrubland<-mcm.habitat[mcm.habitat$Vegetation == "Shrubland",]
write.table(mcm.shrubland[,3:5], file="data.shrubland.csv", sep=",", row.names=F, col.names = F)
# Savanna
mcm.savanna<-mcm.habitat[mcm.habitat$Vegetation == "Savanna",]
write.table(mcm.savanna[,3:5], file="data.savanna.csv", sep=",", row.names=F, col.names = F)
# Tropical seasonal forest
mcm.tropical.seasonal.forest<-mcm.habitat[mcm.habitat$Vegetation == "Tropical seasonal forest",]
write.table(mcm.tropical.seasonal.forest[,3:5], file="data.tropical.seasonal.forest.csv", sep=",", row.names=F, col.names = F)
# Tropical rainforest
mcm.tropical.rainforest<-mcm.habitat[mcm.habitat$Vegetation == "Tropical rainforest",]
write.table(mcm.tropical.rainforest[,3:5], file="data.tropical.rainforest.csv", sep=",", row.names=F, col.names = F)

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
      ## body size disparity
      bs.temp.tropical.rainforest<-abs(data.tropical.rainforest[u.tropical.rainforest,1]-data.tropical.rainforest[v.tropical.rainforest+1,1])
      bs.community.tropical.rainforest<-rbind(bs.community.tropical.rainforest,c(bs.temp.tropical.rainforest))
      bs.tropical.rainforest<-rbind(bs.tropical.rainforest,c(bs.temp.tropical.rainforest))
      ### diet disparity
      dt.temp.tropical.rainforest<-abs(data.tropical.rainforest[u.tropical.rainforest,2]-data.tropical.rainforest[v.tropical.rainforest+1,2])
      dt.community.tropical.rainforest<-rbind(dt.community.tropical.rainforest,c(dt.temp.tropical.rainforest))
      dt.tropical.rainforest<-rbind(dt.tropical.rainforest,c(dt.temp.tropical.rainforest))
      ### locomotor mode disparity
      lm.temp.tropical.rainforest<-abs(data.tropical.rainforest[u.tropical.rainforest,3]-data.tropical.rainforest[v.tropical.rainforest+1,3])
      lm.community.tropical.rainforest<-rbind(lm.community.tropical.rainforest,c(lm.temp.tropical.rainforest))
      lm.tropical.rainforest<-rbind(lm.tropical.rainforest,c(lm.temp.tropical.rainforest))
    }
  }
  mean.bs.community.tropical.rainforest<-rbind(mean.bs.community.tropical.rainforest,c(mean(bs.community.tropical.rainforest[,1])))
  mean.dt.community.tropical.rainforest<-rbind(mean.dt.community.tropical.rainforest,c(mean(dt.community.tropical.rainforest[,1])))
  mean.lm.community.tropical.rainforest<-rbind(mean.lm.community.tropical.rainforest,c(mean(lm.community.tropical.rainforest[,1])))
}

combined.dep.tropical.rainforest<-matrix(c(bs.tropical.rainforest[,1], dt.tropical.rainforest[,1], lm.tropical.rainforest[,1]), ncol=3, byrow = F)

ed.tropical.rainforest.data$Eco.type <- with(ed.tropical.rainforest.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))


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
      ## body size disparity
      bs.temp.tropical.seasonal.forest<-abs(data.tropical.seasonal.forest[u.tropical.seasonal.forest,1]-data.tropical.seasonal.forest[v.tropical.seasonal.forest+1,1])
      bs.community.tropical.seasonal.forest<-rbind(bs.community.tropical.seasonal.forest,c(bs.temp.tropical.seasonal.forest))
      bs.tropical.seasonal.forest<-rbind(bs.tropical.seasonal.forest,c(bs.temp.tropical.seasonal.forest))
      ### diet disparity
      dt.temp.tropical.seasonal.forest<-abs(data.tropical.seasonal.forest[u.tropical.seasonal.forest,2]-data.tropical.seasonal.forest[v.tropical.seasonal.forest+1,2])
      dt.community.tropical.seasonal.forest<-rbind(dt.community.tropical.seasonal.forest,c(dt.temp.tropical.seasonal.forest))
      dt.tropical.seasonal.forest<-rbind(dt.tropical.seasonal.forest,c(dt.temp.tropical.seasonal.forest))
      ### locomotor mode disparity
      lm.temp.tropical.seasonal.forest<-abs(data.tropical.seasonal.forest[u.tropical.seasonal.forest,3]-data.tropical.seasonal.forest[v.tropical.seasonal.forest+1,3])
      lm.community.tropical.seasonal.forest<-rbind(lm.community.tropical.seasonal.forest,c(lm.temp.tropical.seasonal.forest))
      lm.tropical.seasonal.forest<-rbind(lm.tropical.seasonal.forest,c(lm.temp.tropical.seasonal.forest))
    }
  }
  mean.bs.community.tropical.seasonal.forest<-rbind(mean.bs.community.tropical.seasonal.forest,c(mean(bs.community.tropical.seasonal.forest[,1])))
  mean.dt.community.tropical.seasonal.forest<-rbind(mean.dt.community.tropical.seasonal.forest,c(mean(dt.community.tropical.seasonal.forest[,1])))
  mean.lm.community.tropical.seasonal.forest<-rbind(mean.lm.community.tropical.seasonal.forest,c(mean(lm.community.tropical.seasonal.forest[,1])))
}

combined.dep.tropical.seasonal.forest<-matrix(c(bs.tropical.seasonal.forest[,1], dt.tropical.seasonal.forest[,1], lm.tropical.seasonal.forest[,1]), ncol=3, byrow = F)

ed.tropical.seasonal.forest.data$Eco.type <- with(ed.tropical.seasonal.forest.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

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
      ## body size disparity
      bs.temp.savanna<-abs(data.savanna[u.savanna,1]-data.savanna[v.savanna+1,1])
      bs.community.savanna<-rbind(bs.community.savanna,c(bs.temp.savanna))
      bs.savanna<-rbind(bs.savanna,c(bs.temp.savanna))
      ### diet disparity
      dt.temp.savanna<-abs(data.savanna[u.savanna,2]-data.savanna[v.savanna+1,2])
      dt.community.savanna<-rbind(dt.community.savanna,c(dt.temp.savanna))
      dt.savanna<-rbind(dt.savanna,c(dt.temp.savanna))
      ### locomotor mode disparity
      lm.temp.savanna<-abs(data.savanna[u.savanna,3]-data.savanna[v.savanna+1,3])
      lm.community.savanna<-rbind(lm.community.savanna,c(lm.temp.savanna))
      lm.savanna<-rbind(lm.savanna,c(lm.temp.savanna))
    }
  }
  mean.bs.community.savanna<-rbind(mean.bs.community.savanna,c(mean(bs.community.savanna[,1])))
  mean.dt.community.savanna<-rbind(mean.dt.community.savanna,c(mean(dt.community.savanna[,1])))
  mean.lm.community.savanna<-rbind(mean.lm.community.savanna,c(mean(lm.community.savanna[,1])))
}

combined.dep.savanna<-matrix(c(bs.savanna[,1], dt.savanna[,1], lm.savanna[,1]), ncol=3, byrow = F)

ed.savanna.data$Eco.type <- with(ed.savanna.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))


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
      ## body size disparity
      bs.temp.grassland<-abs(data.grassland[u.grassland,1]-data.grassland[v.grassland+1,1])
      bs.community.grassland<-rbind(bs.community.grassland,c(bs.temp.grassland))
      bs.grassland<-rbind(bs.grassland,c(bs.temp.grassland))
      ### diet disparity
      dt.temp.grassland<-abs(data.grassland[u.grassland,2]-data.grassland[v.grassland+1,2])
      dt.community.grassland<-rbind(dt.community.grassland,c(dt.temp.grassland))
      dt.grassland<-rbind(dt.grassland,c(dt.temp.grassland))
      ### locomotor mode disparity
      lm.temp.grassland<-abs(data.grassland[u.grassland,3]-data.grassland[v.grassland+1,3])
      lm.community.grassland<-rbind(lm.community.grassland,c(lm.temp.grassland))
      lm.grassland<-rbind(lm.grassland,c(lm.temp.grassland))
    }
  }
  mean.bs.community.grassland<-rbind(mean.bs.community.grassland,c(mean(bs.community.grassland[,1])))
  mean.dt.community.grassland<-rbind(mean.dt.community.grassland,c(mean(dt.community.grassland[,1])))
  mean.lm.community.grassland<-rbind(mean.lm.community.grassland,c(mean(lm.community.grassland[,1])))
}

combined.dep.grassland<-matrix(c(bs.grassland[,1], dt.grassland[,1], lm.grassland[,1]), ncol=3, byrow = F)

ed.grassland.data$Eco.type <- with(ed.grassland.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

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
      ## body size disparity
      bs.temp.shrubland<-abs(data.shrubland[u.shrubland,1]-data.shrubland[v.shrubland+1,1])
      bs.community.shrubland<-rbind(bs.community.shrubland,c(bs.temp.shrubland))
      bs.shrubland<-rbind(bs.shrubland,c(bs.temp.shrubland))
      ### diet disparity
      dt.temp.shrubland<-abs(data.shrubland[u.shrubland,2]-data.shrubland[v.shrubland+1,2])
      dt.community.shrubland<-rbind(dt.community.shrubland,c(dt.temp.shrubland))
      dt.shrubland<-rbind(dt.shrubland,c(dt.temp.shrubland))
      ### locomotor mode disparity
      lm.temp.shrubland<-abs(data.shrubland[u.shrubland,3]-data.shrubland[v.shrubland+1,3])
      lm.community.shrubland<-rbind(lm.community.shrubland,c(lm.temp.shrubland))
      lm.shrubland<-rbind(lm.shrubland,c(lm.temp.shrubland))
    }
  }
  mean.bs.community.shrubland<-rbind(mean.bs.community.shrubland,c(mean(bs.community.shrubland[,1])))
  mean.dt.community.shrubland<-rbind(mean.dt.community.shrubland,c(mean(dt.community.shrubland[,1])))
  mean.lm.community.shrubland<-rbind(mean.lm.community.shrubland,c(mean(lm.community.shrubland[,1])))
}

combined.dep.shrubland<-matrix(c(bs.shrubland[,1], dt.shrubland[,1], lm.shrubland[,1]), ncol=3, byrow = F)

ed.shrubland.data$Eco.type <- with(ed.shrubland.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))

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
      ## body size disparity
      bs.temp.desert<-abs(data.desert[u.desert,1]-data.desert[v.desert+1,1])
      bs.community.desert<-rbind(bs.community.desert,c(bs.temp.desert))
      bs.desert<-rbind(bs.desert,c(bs.temp.desert))
      ### diet disparity
      dt.temp.desert<-abs(data.desert[u.desert,2]-data.desert[v.desert+1,2])
      dt.community.desert<-rbind(dt.community.desert,c(dt.temp.desert))
      dt.desert<-rbind(dt.desert,c(dt.temp.desert))
      ### locomotor mode disparity
      lm.temp.desert<-abs(data.desert[u.desert,3]-data.desert[v.desert+1,3])
      lm.community.desert<-rbind(lm.community.desert,c(lm.temp.desert))
      lm.desert<-rbind(lm.desert,c(lm.temp.desert))
    }
  }
  mean.bs.community.desert<-rbind(mean.bs.community.desert,c(mean(bs.community.desert[,1])))
  mean.dt.community.desert<-rbind(mean.dt.community.desert,c(mean(dt.community.desert[,1])))
  mean.lm.community.desert<-rbind(mean.lm.community.desert,c(mean(lm.community.desert[,1])))
}

combined.dep.desert<-matrix(c(bs.desert[,1], dt.desert[,1], lm.desert[,1]), ncol=3, byrow = F)

ed.desert.data$Eco.type <- with(ed.desert.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))


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
      ## body size disparity
      bs.temp.temperate.forest<-abs(data.temperate.forest[u.temperate.forest,1]-data.temperate.forest[v.temperate.forest+1,1])
      bs.community.temperate.forest<-rbind(bs.community.temperate.forest,c(bs.temp.temperate.forest))
      bs.temperate.forest<-rbind(bs.temperate.forest,c(bs.temp.temperate.forest))
      ### diet disparity
      dt.temp.temperate.forest<-abs(data.temperate.forest[u.temperate.forest,2]-data.temperate.forest[v.temperate.forest+1,2])
      dt.community.temperate.forest<-rbind(dt.community.temperate.forest,c(dt.temp.temperate.forest))
      dt.temperate.forest<-rbind(dt.temperate.forest,c(dt.temp.temperate.forest))
      ### locomotor mode disparity
      lm.temp.temperate.forest<-abs(data.temperate.forest[u.temperate.forest,3]-data.temperate.forest[v.temperate.forest+1,3])
      lm.community.temperate.forest<-rbind(lm.community.temperate.forest,c(lm.temp.temperate.forest))
      lm.temperate.forest<-rbind(lm.temperate.forest,c(lm.temp.temperate.forest))
    }
  }
  mean.bs.community.temperate.forest<-rbind(mean.bs.community.temperate.forest,c(mean(bs.community.temperate.forest[,1])))
  mean.dt.community.temperate.forest<-rbind(mean.dt.community.temperate.forest,c(mean(dt.community.temperate.forest[,1])))
  mean.lm.community.temperate.forest<-rbind(mean.lm.community.temperate.forest,c(mean(lm.community.temperate.forest[,1])))
}

combined.dep.temperate.forest<-matrix(c(bs.temperate.forest[,1], dt.temperate.forest[,1], lm.temperate.forest[,1]), ncol=3, byrow = F)

ed.temperate.forest.data$Eco.type <- with(ed.temperate.forest.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))


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
      ## body size disparity
      bs.temp.boreal<-abs(data.boreal[u.boreal,1]-data.boreal[v.boreal+1,1])
      bs.community.boreal<-rbind(bs.community.boreal,c(bs.temp.boreal))
      bs.boreal<-rbind(bs.boreal,c(bs.temp.boreal))
      ### diet disparity
      dt.temp.boreal<-abs(data.boreal[u.boreal,2]-data.boreal[v.boreal+1,2])
      dt.community.boreal<-rbind(dt.community.boreal,c(dt.temp.boreal))
      dt.boreal<-rbind(dt.boreal,c(dt.temp.boreal))
      ### locomotor mode disparity
      lm.temp.boreal<-abs(data.boreal[u.boreal,3]-data.boreal[v.boreal+1,3])
      lm.community.boreal<-rbind(lm.community.boreal,c(lm.temp.boreal))
      lm.boreal<-rbind(lm.boreal,c(lm.temp.boreal))
    }
  }
  mean.bs.community.boreal<-rbind(mean.bs.community.boreal,c(mean(bs.community.boreal[,1])))
  mean.dt.community.boreal<-rbind(mean.dt.community.boreal,c(mean(dt.community.boreal[,1])))
  mean.lm.community.boreal<-rbind(mean.lm.community.boreal,c(mean(lm.community.boreal[,1])))
}

combined.dep.boreal<-matrix(c(bs.boreal[,1], dt.boreal[,1], lm.boreal[,1]), ncol=3, byrow = F)

##
ed.boreal.data$Eco.type <- with(ed.boreal.data, paste0(BS.Rank, Diet.Rank, Locomotor.Rank))


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

write.csv(vegetation.dep, file = "EDisp_vegetation_statistic_24July2017.csv")
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

write.csv(vegetation.ed, file = "ERich_vegetation_statistic_24July2017.csv")

### Plots of proportoin data of eight vegetation types
par(mfcol=c(8,3),
    oma=c(1,1,0,1), 
    mar=c(2,2,1,0), 
    mgp=c(1,0.6,0),
    tck=-0.02)

# body size
barplot(table(ed.tropical.rainforest.data[,1])/sum(table(ed.tropical.rainforest.data[,1])), ylim=c(0,0.6), space=1, col="red") # tropical rainforest
barplot(table(ed.tropical.seasonal.forest.data[,1])/sum(table(ed.tropical.seasonal.forest.data[,1])), ylim=c(0,0.6), space=1, col="orangered") # tropical seasonal foest
barplot(table(ed.savanna.data[,1])/sum(table(ed.savanna.data[,1])), ylim=c(0,0.6), space=1, col="firebrick1") # savanna
barplot(table(ed.grassland.data[,1])/sum(table(ed.grassland.data[,1])), ylim=c(0,0.6), space=1, col="goldenrod") # grassland
barplot(table(ed.shrubland.data[,1])/sum(table(ed.shrubland.data[,1])), ylim=c(0,0.6), space=1, col="darkgoldenrod1") # shrubland
barplot(table(ed.desert.data[,1])/sum(table(ed.desert.data[,1])), ylim=c(0,0.6), space=1, col="khaki") # desert
barplot(table(ed.temperate.forest.data[,1])/sum(table(ed.temperate.forest.data[,1])), ylim=c(0,0.6), space=1, col="green3") # temperate forest
barplot(table(ed.boreal.data[,1])/sum(table(ed.boreal.data[,1])), ylim=c(0,0.6), space=1, col="darkviolet") # boreal forest


# diet
barplot(table(ed.tropical.rainforest.data[,2])/sum(table(ed.tropical.rainforest.data[,2])), ylim=c(0,0.6), space=1, col="red") # tropical rainforest
barplot(table(ed.tropical.seasonal.forest.data[,2])/sum(table(ed.tropical.seasonal.forest.data[,2])), ylim=c(0,0.6), space=1, col="orangered") # tropical seasonal foest
barplot(table(ed.savanna.data[,2])/sum(table(ed.savanna.data[,2])), ylim=c(0,0.6), space=1, col="firebrick1") # savanna
barplot(table(ed.grassland.data[,2])/sum(table(ed.grassland.data[,2])), ylim=c(0,0.6), space=1, col="goldenrod") # grassland
barplot(table(ed.shrubland.data[,2])/sum(table(ed.shrubland.data[,2])), ylim=c(0,0.6), space=1, col="darkgoldenrod1") # shrubland
barplot(table(ed.desert.data[,2])/sum(table(ed.desert.data[,2])), ylim=c(0,0.6), space=1, col="khaki") # desert
barplot(table(ed.temperate.forest.data[,2])/sum(table(ed.temperate.forest.data[,2])), ylim=c(0,0.6), space=1, col="green3") # temperate forest
barplot(table(ed.boreal.data[,2])/sum(table(ed.boreal.data[,2])), ylim=c(0,0.6), space=1, col="darkviolet") # boreal forest

# locomotor mode
barplot(table(ed.tropical.rainforest.data[,3])/sum(table(ed.tropical.rainforest.data[,3])), ylim=c(0,0.6), space=1, col="red") # tropical rainforest
barplot(table(ed.tropical.seasonal.forest.data[,3])/sum(table(ed.tropical.seasonal.forest.data[,3])), ylim=c(0,0.6), space=1, col="orangered") # tropical seasonal foest
barplot(table(ed.savanna.data[,3])/sum(table(ed.savanna.data[,3])), ylim=c(0,0.6), space=1, col="firebrick1") # savanna
barplot(table(ed.grassland.data[,3])/sum(table(ed.grassland.data[,3])), ylim=c(0,0.6), space=1, col="goldenrod") # grassland
barplot(table(ed.shrubland.data[,3])/sum(table(ed.shrubland.data[,3])), ylim=c(0,0.6), space=1, col="darkgoldenrod1") # shrubland
barplot(table(ed.desert.data[,3])/sum(table(ed.desert.data[,3])), ylim=c(0,0.6), space=1, col="khaki") # desert
barplot(table(ed.temperate.forest.data[,3])/sum(table(ed.temperate.forest.data[,3])), ylim=c(0,0.6), space=1, col="green3") # temperate forest
barplot(table(ed.boreal.data[,3])/sum(table(ed.boreal.data[,3])), ylim=c(0,0.6), space=1, col="darkviolet") # boreal forest

###########################################
### Combine all proporiton data together
############################################

# body size
bs.proportion.all <- rbind(table(ed.close.data[,1])/sum(table(ed.close.data[,1])),
      table(ed.open.data[,1])/sum(table(ed.open.data[,1])),
      table(ed.tropical.climate.data[,1])/sum(table(ed.tropical.climate.data[,1])), 
      table(ed.arid.climate.data[,1])/sum(table(ed.arid.climate.data[,1])), 
      table(ed.temperate.climate.data[,1])/sum(table(ed.temperate.climate.data[,1])),
      table(ed.cold.climate.data[,1])/sum(table(ed.cold.climate.data[,1])), 
      table(ed.tropical.rainforest.data[,1])/sum(table(ed.tropical.rainforest.data[,1])),
      table(ed.tropical.seasonal.forest.data[,1])/sum(table(ed.tropical.seasonal.forest.data[,1])),
      table(ed.savanna.data[,1])/sum(table(ed.savanna.data[,1])), 
      table(ed.grassland.data[,1])/sum(table(ed.grassland.data[,1])),
      table(ed.shrubland.data[,1])/sum(table(ed.shrubland.data[,1])),
      table(ed.desert.data[,1])/sum(table(ed.desert.data[,1])),
      table(ed.temperate.forest.data[,1])/sum(table(ed.temperate.forest.data[,1])),
      table(ed.boreal.data[,1])/sum(table(ed.boreal.data[,1])))
# diet
dt.proportion.all <-rbind(table(ed.close.data[,2])/sum(table(ed.close.data[,2])), 
      table(ed.open.data[,2])/sum(table(ed.open.data[,2])),
      table(ed.tropical.climate.data[,2])/sum(table(ed.tropical.climate.data[,2])),
      table(ed.arid.climate.data[,2])/sum(table(ed.arid.climate.data[,2])),
      table(ed.temperate.climate.data[,2])/sum(table(ed.temperate.climate.data[,2])),
      table(ed.cold.climate.data[,2])/sum(table(ed.cold.climate.data[,2])),
      table(ed.tropical.rainforest.data[,2])/sum(table(ed.tropical.rainforest.data[,2])), 
      table(ed.tropical.seasonal.forest.data[,2])/sum(table(ed.tropical.seasonal.forest.data[,2])), 
      table(ed.savanna.data[,2])/sum(table(ed.savanna.data[,2])),
      table(ed.grassland.data[,2])/sum(table(ed.grassland.data[,2])),
      table(ed.shrubland.data[,2])/sum(table(ed.shrubland.data[,2])), 
      table(ed.desert.data[,2])/sum(table(ed.desert.data[,2])),
      table(ed.temperate.forest.data[,2])/sum(table(ed.temperate.forest.data[,2])),
      table(ed.boreal.data[,2])/sum(table(ed.boreal.data[,2])))

# locomotor mode
lm.proportion.all <-rbind(table(ed.close.data[,3])/sum(table(ed.close.data[,3])), 
      table(ed.open.data[,3])/sum(table(ed.open.data[,3])),
      table(ed.tropical.climate.data[,3])/sum(table(ed.tropical.climate.data[,3])),
      table(ed.arid.climate.data[,3])/sum(table(ed.arid.climate.data[,3])),
      table(ed.temperate.climate.data[,3])/sum(table(ed.temperate.climate.data[,3])),
      table(ed.cold.climate.data[,3])/sum(table(ed.cold.climate.data[,3])),
      table(ed.tropical.rainforest.data[,3])/sum(table(ed.tropical.rainforest.data[,3])), 
      table(ed.tropical.seasonal.forest.data[,3])/sum(table(ed.tropical.seasonal.forest.data[,3])), 
      table(ed.savanna.data[,3])/sum(table(ed.savanna.data[,3])),
      table(ed.grassland.data[,3])/sum(table(ed.grassland.data[,3])),
      table(ed.shrubland.data[,3])/sum(table(ed.shrubland.data[,3])),
      table(ed.desert.data[,3])/sum(table(ed.desert.data[,3])),
      table(ed.temperate.forest.data[,3])/sum(table(ed.temperate.forest.data[,3])),
      table(ed.boreal.data[,3])/sum(table(ed.boreal.data[,3])))

proportion.all <- cbind(bs.proportion.all, dt.proportion.all, lm.proportion.all)
proportion.all <- proportion.all*100

colnames(proportion.all) <- c("BS1", "BS2", "BS3", "BS4", "BS5", 
                              "DT1", "DT2", "DT3", "DT4", "DT5", "DT6",
                              "LM1", "LM2", "LM3", "LM4", "LM5", "LM6", "LM7", "LM8")
rownames(proportion.all) <- c("close", "open",
                              "tropical", "arid", "temperate", "cold",
                              "tropical rain forest", "tropical seasonal forest", "savanna", "grassland", 
                              "shrubland", "desert", "temperate forest", "boreal forest")
########### this data need to clean up because it missing some data that not reserved during the table function
# write.csv(proportion.all, file = "functional trait percentage_10August2017.csv", sep = ",")  

###########################################################
#### Edisp contribution of each eoclogical paratmer
###########################################################

par(mfcol=c(1,1))

#### stack plots
tropical.rainforest.type<-rbind(sum(bs.tropical.rainforest[,1])/sum(combined.dep.tropical.rainforest), sum(dt.tropical.rainforest[,1])/sum(combined.dep.tropical.rainforest), 
                     sum(lm.tropical.rainforest[,1])/sum(combined.dep.tropical.rainforest))
tropical.seasonal.forest.type<-rbind(sum(bs.tropical.seasonal.forest[,1])/sum(combined.dep.tropical.seasonal.forest), sum(dt.tropical.seasonal.forest[,1])/sum(combined.dep.tropical.seasonal.forest), 
                                sum(lm.tropical.seasonal.forest[,1])/sum(combined.dep.tropical.seasonal.forest))
desert.type<-rbind(sum(bs.desert[,1])/sum(combined.dep.desert), sum(dt.desert[,1])/sum(combined.dep.desert), 
                               sum(lm.desert[,1])/sum(combined.dep.desert))
savanna.type<-rbind(sum(bs.savanna[,1])/sum(combined.dep.savanna), sum(dt.savanna[,1])/sum(combined.dep.savanna), 
                     sum(lm.savanna[,1])/sum(combined.dep.savanna))
grassland.type<-rbind(sum(bs.grassland[,1])/sum(combined.dep.grassland), sum(dt.grassland[,1])/sum(combined.dep.grassland), 
                     sum(lm.grassland[,1])/sum(combined.dep.grassland))
shrubland.type<-rbind(sum(bs.shrubland[,1])/sum(combined.dep.shrubland), sum(dt.shrubland[,1])/sum(combined.dep.shrubland), 
                     sum(lm.shrubland[,1])/sum(combined.dep.shrubland))
temperate.forest.type<-rbind(sum(bs.temperate.forest[,1])/sum(combined.dep.temperate.forest), sum(dt.temperate.forest[,1])/sum(combined.dep.temperate.forest), 
                     sum(lm.temperate.forest[,1])/sum(combined.dep.temperate.forest))
boreal.type<-rbind(sum(bs.boreal[,1])/sum(combined.dep.boreal), sum(dt.boreal[,1])/sum(combined.dep.boreal), 
                     sum(lm.boreal[,1])/sum(combined.dep.boreal))


#### stack plots
combined.vegetation.type<-cbind(tropical.rainforest.type*mean(dep.tropical.rainforest[,1]),
                                tropical.seasonal.forest.type*mean(dep.tropical.seasonal.forest[,1]),
                                savanna.type*mean(dep.savanna[,1]),
                                grassland.type*mean(dep.grassland[,1]),
                                shrubland.type*mean(dep.shrubland[,1]),
                                temperate.forest.type*mean(dep.temperate.forest[,1]),
                                desert.type*mean(dep.desert[,1]),
                                boreal.type*mean(dep.boreal[,1]))
barplot(combined.vegetation.type, space=1, col=c("dodgerblue3","green3","darkorange"))

### Habitat, climatic, and vegetation types together
par(mfrow=c(1,1),oma=c(3,3,3,1),mar=c(1,1,2,1),
    mgp=c(3,0.5,0),tck=-0.005,cex=0.8,cex.axis=0.7)

sd.all<-cbind(sd(dep.close[,1]),
              sd(dep.open[,1]),
              sd(dep.tropical.climate[,1]),
              sd(dep.arid.climate[,1]),
              sd(dep.temperate.climate[,1]),
              sd(dep.cold.climate[,1]),
              sd(dep.tropical.rainforest[,1]),
              sd(dep.tropical.seasonal.forest[,1]),
              sd(dep.savanna[,1]),
              sd(dep.shrubland[,1]),
              sd(dep.desert[,1]),
              sd(dep.grassland[,1]),
              sd(dep.temperate.forest[,1]),
              sd(dep.boreal[,1]))

length.all<-t(matrix(cbind(length(dep.close[,1]),
                           length(dep.open[,1]),
                           length(dep.tropical.climate[,1]),
                           length(dep.arid.climate[,1]),
                           length(dep.temperate.climate[,1]),
                           length(dep.cold.climate[,1]),
                           length(dep.tropical.rainforest[,1]),
                           length(dep.tropical.seasonal.forest[,1]),
                           length(dep.savanna[,1]),
                           length(dep.shrubland[,1]),
                           length(dep.desert[,1]),
                           length(dep.grassland[,1]),
                           length(dep.temperate.forest[,1]),
                           length(dep.boreal[,1]))))

combined.all<-cbind(mean(dep.close[,1]),
                    mean(dep.open[,1]),
                    mean(dep.tropical.climate[,1]),
                    mean(dep.arid.climate[,1]),
                    mean(dep.temperate.climate[,1]),
                    mean(dep.cold.climate[,1]),
                    mean(dep.tropical.rainforest[,1]),
                    mean(dep.tropical.seasonal.forest[,1]),
                    mean(dep.savanna[,1]),
                    mean(dep.shrubland[,1]),
                    mean(dep.desert[,1]),
                    mean(dep.grassland[,1]),
                    mean(dep.temperate.forest[,1]),
                    mean(dep.boreal[,1]))

combined.all.type<-cbind(close.type*mean(dep.close[,1]),
                         open.type*mean(dep.open[,1]),
                         tropical.climatic.type*mean(dep.tropical.climate[,1]),
                         arid.climatic.type*mean(dep.arid.climate[,1]),
                         temperate.climatic.type*mean(dep.temperate.climate[,1]),
                         cold.climatic.type*mean(dep.cold.climate[,1]),
                         tropical.rainforest.type*mean(dep.tropical.rainforest[,1]),
                         tropical.seasonal.forest.type*mean(dep.tropical.seasonal.forest[,1]),
                         savanna.type*mean(dep.savanna[,1]),
                         shrubland.type*mean(dep.shrubland[,1]),
                         desert.type*mean(dep.desert[,1]),
                         grassland.type*mean(dep.grassland[,1]),
                         temperate.forest.type*mean(dep.temperate.forest[,1]),
                         boreal.type*mean(dep.boreal[,1]))

barx<-barplot(combined.all.type, 
              space=1, col=c("burlywood4","cyan4","gold1"), 
              ylim=c(0, 5),
              names.arg=c("close","open","tropical","arid","temperate","cold","tropical rainforest",
                    "tropical seasonal forest","savanna","shrubland",
                    "temperate desert","grassland","temperate forest","boreal forest"))

error.bar(barx, combined.all, 1.96*sd.all/sqrt(length.all))

################################################
## Dissimilarity among different vegetations
################################################

## tropical rain forest agaisnt others
dis.tropicalrainforest.tropicalseasaonlforest<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                   unique(ed.tropical.seasonal.forest.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                   unique(ed.tropical.seasonal.forest.data$Eco.type)))
dis.tropicalrainforest.savanna<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                   unique(ed.savanna.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                   unique(ed.savanna.data$Eco.type)))
dis.tropicalrainforest.grassland<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                 unique(ed.grassland.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                 unique(ed.grassland.data$Eco.type)))
dis.tropicalrainforest.shrubland<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                 unique(ed.shrubland.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                 unique(ed.shrubland.data$Eco.type)))
dis.tropicalrainforest.desert<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                       unique(ed.desert.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                       unique(ed.desert.data$Eco.type)))
dis.tropicalrainforest.temperateforest<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                    unique(ed.temperate.forest.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                    unique(ed.temperate.forest.data$Eco.type)))
dis.tropicalrainforest.borealforest<-1-length(intersect(unique(ed.tropical.rainforest.data$Eco.type),
                                                 unique(ed.boreal.data$Eco.type)))/length(union(unique(ed.tropical.rainforest.data$Eco.type),
                                                                                                 unique(ed.boreal.data$Eco.type)))
## torpical seasoanl forest against others
dis.tropicalseasonalforest.savanna<-1-length(intersect(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                   unique(ed.savanna.data$Eco.type)))/length(union(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                                                                   unique(ed.savanna.data$Eco.type)))
dis.tropicalseasonalforest.grassland<-1-length(intersect(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                     unique(ed.grassland.data$Eco.type)))/length(union(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                                                                       unique(ed.grassland.data$Eco.type)))
dis.tropicalseasonalforest.shrubland<-1-length(intersect(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                     unique(ed.shrubland.data$Eco.type)))/length(union(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                                                                       unique(ed.shrubland.data$Eco.type)))
dis.tropicalseasonalforest.desert<-1-length(intersect(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                  unique(ed.desert.data$Eco.type)))/length(union(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                                                                 unique(ed.desert.data$Eco.type)))
dis.tropicalseasonalforest.temperateforest<-1-length(intersect(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                           unique(ed.temperate.forest.data$Eco.type)))/length(union(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                                                                                    unique(ed.temperate.forest.data$Eco.type)))
dis.tropicalseasonalforest.borealforest<-1-length(intersect(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                          unique(ed.boreal.data$Eco.type)))/length(union(unique(ed.tropical.seasonal.forest.data$Eco.type),
                                                                                                         unique(ed.boreal.data$Eco.type)))
## savanna against others
dis.savanna.grassland<-1-length(intersect(unique(ed.savanna.data$Eco.type),
                                          unique(ed.grassland.data$Eco.type)))/length(union(unique(ed.savanna.data$Eco.type),
                                                                                          unique(ed.grassland.data$Eco.type)))
dis.savanna.shrubland<-1-length(intersect(unique(ed.savanna.data$Eco.type),
                                          unique(ed.shrubland.data$Eco.type)))/length(union(unique(ed.savanna.data$Eco.type),
                                                                                            unique(ed.shrubland.data$Eco.type)))
dis.savanna.desert<-1-length(intersect(unique(ed.savanna.data$Eco.type),
                                          unique(ed.desert.data$Eco.type)))/length(union(unique(ed.savanna.data$Eco.type),
                                                                                            unique(ed.desert.data$Eco.type)))
dis.savanna.temperateforest<-1-length(intersect(unique(ed.savanna.data$Eco.type),
                                          unique(ed.temperate.forest.data$Eco.type)))/length(union(unique(ed.savanna.data$Eco.type),
                                                                                            unique(ed.temperate.forest.data$Eco.type)))
dis.savanna.borealforest<-1-length(intersect(unique(ed.savanna.data$Eco.type),
                                          unique(ed.boreal.data$Eco.type)))/length(union(unique(ed.savanna.data$Eco.type),
                                                                                            unique(ed.boreal.data$Eco.type)))

## grassland against others
                                                       
dis.grassland.shrubland<-1-length(intersect(unique(ed.grassland.data$Eco.type),
                                            unique(ed.shrubland.data$Eco.type)))/length(union(unique(ed.grassland.data$Eco.type),
                                                                                              unique(ed.shrubland.data$Eco.type)))
dis.grassland.desert<-1-length(intersect(unique(ed.grassland.data$Eco.type),
                                            unique(ed.desert.data$Eco.type)))/length(union(unique(ed.grassland.data$Eco.type),
                                                                                              unique(ed.desert.data$Eco.type)))
dis.grassland.temperateforest<-1-length(intersect(unique(ed.grassland.data$Eco.type),
                                            unique(ed.temperate.forest.data$Eco.type)))/length(union(unique(ed.grassland.data$Eco.type),
                                                                                              unique(ed.temperate.forest.data$Eco.type)))
dis.grassland.borealforest<-1-length(intersect(unique(ed.grassland.data$Eco.type),
                                            unique(ed.boreal.data$Eco.type)))/length(union(unique(ed.grassland.data$Eco.type),
                                                                                              unique(ed.boreal.data$Eco.type)))
## shrubland against others
dis.shrubland.desert<-1-length(intersect(unique(ed.shrubland.data$Eco.type),
                                            unique(ed.desert.data$Eco.type)))/length(union(unique(ed.shrubland.data$Eco.type),
                                                                                              unique(ed.desert.data$Eco.type)))
dis.shrubland.temperateforest<-1-length(intersect(unique(ed.shrubland.data$Eco.type),
                                         unique(ed.temperate.forest.data$Eco.type)))/length(union(unique(ed.shrubland.data$Eco.type),
                                                                                        unique(ed.temperate.forest.data$Eco.type)))
dis.shrubland.borealfoest<-1-length(intersect(unique(ed.shrubland.data$Eco.type),
                                         unique(ed.boreal.data$Eco.type)))/length(union(unique(ed.shrubland.data$Eco.type),
                                                                                        unique(ed.boreal.data$Eco.type)))

## desert against others
dis.desert.temperateforest<-1-length(intersect(unique(ed.temperate.forest.data$Eco.type),
                                         unique(ed.desert.data$Eco.type)))/length(union(unique(ed.temperate.forest.data$Eco.type),
                                                                                        unique(ed.desert.data$Eco.type)))
dis.desert.borealforest<-1-length(intersect(unique(ed.boreal.data$Eco.type),
                                         unique(ed.desert.data$Eco.type)))/length(union(unique(ed.boreal.data$Eco.type),
                                                                                        unique(ed.desert.data$Eco.type)))
## temperate forest against other
dis.temperateforest.borealforest<-1-length(intersect(unique(ed.temperate.forest.data$Eco.type),
                                               unique(ed.boreal.data$Eco.type)))/length(union(unique(ed.temperate.forest.data$Eco.type),
                                                                                              unique(ed.boreal.data$Eco.type)))
dis.tropicalrainforest.tropicalseasaonlforest
dis.tropicalrainforest.savanna
dis.tropicalrainforest.grassland
dis.tropicalrainforest.shrubland
dis.tropicalrainforest.desert
dis.tropicalrainforest.temperateforest
dis.tropicalrainforest.borealforest
## torpical seasoanl forest against others
dis.tropicalseasonalforest.savanna
dis.tropicalseasonalforest.grassland
dis.tropicalseasonalforest.shrubland
dis.tropicalseasonalforest.desert
dis.tropicalseasonalforest.temperateforest
dis.tropicalseasonalforest.borealforest
## savanna against others
dis.savanna.grassland
dis.savanna.shrubland
dis.savanna.desert
dis.savanna.temperateforest
dis.savanna.borealforest
## grassland against others
dis.grassland.shrubland
dis.grassland.desert
dis.grassland.temperateforest
dis.grassland.borealforest
## shrubland against others
dis.shrubland.desert
dis.shrubland.temperateforest
dis.shrubland.borealfoest
## desert against others
dis.desert.temperateforest
dis.desert.borealforest
## temperate forest against others
dis.temperateforest.borealforest

###################################
## Ecological Richness (ERich) of all environmental types
###################################

par(mfrow=c(1,1),oma=c(3,3,1,1),mar=c(1,1,1,1),
    mgp=c(3,0.5,0),tck=-0.005,cex=0.8,cex.axis=0.7)

sd.ed.all<-cbind(sd(ed.close[,1]),
              sd(ed.open[,1]),
              sd(ed.tropical.climate[,1]),
              sd(ed.arid.climate[,1]),
              sd(ed.temperate.climate[,1]),
              sd(ed.cold.climate[,1]),
              sd(ed.tropical.rainforest[,1]),
              sd(ed.tropical.seasonal.forest[,1]),
              sd(ed.savanna[,1]),
              sd(ed.shrubland[,1]),
              sd(ed.desert[,1]),
              sd(ed.grassland[,1]),
              sd(ed.temperate.forest[,1]),
              sd(ed.boreal[,1]))

length.ed.all<-t(matrix(cbind(length(ed.close[,1]),
                           length(ed.open[,1]),
                           length(ed.tropical.climate[,1]),
                           length(ed.arid.climate[,1]),
                           length(ed.temperate.climate[,1]),
                           length(ed.cold.climate[,1]),
                           length(ed.tropical.rainforest[,1]),
                           length(ed.tropical.seasonal.forest[,1]),
                           length(ed.savanna[,1]),
                           length(ed.shrubland[,1]),
                           length(ed.desert[,1]),
                           length(ed.grassland[,1]),
                           length(ed.temperate.forest[,1]),
                           length(ed.boreal[,1]))))

combined.ed.all<-cbind(mean(ed.close[,1]),
                    mean(ed.open[,1]),
                    mean(ed.tropical.climate[,1]),
                    mean(ed.arid.climate[,1]),
                    mean(ed.temperate.climate[,1]),
                    mean(ed.cold.climate[,1]),
                    mean(ed.tropical.rainforest[,1]),
                    mean(ed.tropical.seasonal.forest[,1]),
                    mean(ed.savanna[,1]),
                    mean(ed.shrubland[,1]),
                    mean(ed.desert[,1]),
                    mean(ed.grassland[,1]),
                    mean(ed.temperate.forest[,1]),
                    mean(ed.boreal[,1]))

## all together
barx <- barplot(combined.ed.all, 
                space=1, col=c("grey70"), 
                ylim=c(0, 15),
                names.arg=c("close","open","tropical","arid","temperate","cold","tropical rainforest",
                            "tropical seasonal forest","savanna","shrubland",
                            "temperate desert","grassland","temperate forest","boreal forest"),
                las=1, cex.names = 0.8, axis.lty = 1)
  
## plot the error bar
error.bar(barx, combined.ed.all, 1.96*sd.ed.all/sqrt(length.ed.all))
mtext(outer=F, side=2, line=1.5, text="Mean Diversity",cex=0.8)


#####################################################
#### Pair-wise t-test between different environments
#####################################################

## EDisp
## habitats
dep.t.test.close_open <- t.test(dep.close[,1], dep.open[,1])
## climates
dep.t.test.tropical_arid <- t.test(dep.tropical.climate[,1], dep.arid.climate[,1])
dep.t.test.tropical_temperate <- t.test(dep.tropical.climate[,1], dep.temperate.climate[,1])
dep.t.test.tropical_cold <- t.test(dep.tropical.climate[,1], dep.cold.climate[,1])
dep.t.test.arid_temperate <- t.test(dep.arid.climate[,1], dep.temperate.climate[,1])
dep.t.test.arid_cold <- t.test(dep.arid.climate[,1], dep.cold.climate[,1])
dep.t.test.temperate_cold <- t.test(dep.temperate.climate[,1], dep.cold.climate[,1])
### vegetations
dep.t.test.tropical.rainforest_tropical.seasonalforest <- t.test(dep.tropical.rainforest[,1], dep.tropical.seasonal.forest[,1])
dep.t.test.tropical.rainforest_savanna<-t.test(dep.tropical.rainforest[,1], dep.savanna[,1])
dep.t.test.tropical.rainforest_grassland <- t.test(dep.tropical.rainforest[,1], dep.grassland[,1])
dep.t.test.tropical.rainforest_shrubland <- t.test(dep.tropical.rainforest[,1], dep.shrubland[,1])
dep.t.test.tropical.rainforest_desert <- t.test(dep.tropical.rainforest[,1], dep.desert[,1])
dep.t.test.tropical.rainforest_temperate.forest <- t.test(dep.tropical.rainforest[,1], dep.temperate.forest[,1])
dep.t.test.tropical.rainforest_boreal <- t.test(dep.tropical.rainforest[,1], dep.boreal[,1])
dep.t.test.tropical.seasonalforest_savanna<- t.test(dep.tropical.seasonal.forest[,1], dep.savanna[,1])
dep.t.test.tropical.seasonalforest_grassland <- t.test(dep.tropical.seasonal.forest[,1], dep.grassland[,1])
dep.t.test.tropical.seasonalforest_shrubland <- t.test(dep.tropical.seasonal.forest[,1], dep.shrubland[,1])
dep.t.test.tropical.seasonalforest_desert <- t.test(dep.tropical.seasonal.forest[,1], dep.desert[,1])
dep.t.test.tropical.seasonalforest_temperate.forest <- t.test(dep.tropical.seasonal.forest[,1], dep.temperate.forest[,1])
dep.t.test.tropical.seasonalforest_boreal <- t.test(dep.tropical.seasonal.forest[,1], dep.boreal[,1])
dep.t.test.savanna_grassland <- t.test(dep.savanna[,1], dep.grassland[,1])
dep.t.test.savanna_shrubland <- t.test(dep.savanna[,1], dep.shrubland[,1])
dep.t.test.savanna_desert <- t.test(dep.savanna[,1], dep.desert[,1])
dep.t.test.savanna_temperate.forest <- t.test(dep.savanna[,1], dep.temperate.forest[,1])
dep.t.test.savanna_boreal <- t.test(dep.savanna[,1], dep.boreal[,1])
dep.t.test.grassland_shrubland <- t.test(dep.grassland[,1], dep.shrubland[,1])
dep.t.test.grassland_desert <- t.test(dep.grassland[,1], dep.desert[,1])
dep.t.test.grassland_temperate.forest <- t.test(dep.grassland[,1], dep.temperate.forest[,1])
dep.t.test.grassland_boreal <- t.test(dep.grassland[,1], dep.boreal[,1])
dep.t.test.shrubland_desert <- t.test(dep.shrubland[,1], dep.desert[,1])
dep.t.test.shrubland_temperate.forest <- t.test(dep.shrubland[,1], dep.temperate.forest[,1])
dep.t.test.shrubland_boreal <- t.test(dep.shrubland[,1], dep.boreal[,1])
dep.t.test.desert_temperate.forest <- t.test(dep.desert[,1], dep.temperate.forest[,1])
dep.t.test.desert_boreal <- t.test(dep.desert[,1], dep.boreal[,1])
dep.t.test.temperate_boreal <- t.test(dep.temperate.forest[,1], dep.boreal[,1])

pairwise.t.test.EDisp <- rbind(dep.t.test.close_open,
                               dep.t.test.tropical_arid,
                               dep.t.test.tropical_temperate,
                               dep.t.test.tropical_cold,
                               dep.t.test.arid_temperate,
                               dep.t.test.arid_cold,
                               dep.t.test.temperate_cold,
                               dep.t.test.tropical.rainforest_tropical.seasonalforest,
                               dep.t.test.tropical.rainforest_savanna,
                               dep.t.test.tropical.rainforest_grassland,
                               dep.t.test.tropical.rainforest_shrubland,
                               dep.t.test.tropical.rainforest_desert,
                               dep.t.test.tropical.rainforest_temperate.forest,
                               dep.t.test.tropical.rainforest_boreal,
                               dep.t.test.tropical.seasonalforest_savanna,
                               dep.t.test.tropical.seasonalforest_grassland,
                               dep.t.test.tropical.seasonalforest_shrubland,
                               dep.t.test.tropical.seasonalforest_desert,
                               dep.t.test.tropical.seasonalforest_temperate.forest,
                               dep.t.test.tropical.seasonalforest_boreal,
                               dep.t.test.savanna_grassland,
                               dep.t.test.savanna_shrubland,
                               dep.t.test.savanna_desert,
                               dep.t.test.savanna_temperate.forest,
                               dep.t.test.savanna_boreal,
                               dep.t.test.grassland_shrubland,
                               dep.t.test.grassland_desert,
                               dep.t.test.grassland_temperate.forest,
                               dep.t.test.grassland_boreal,
                               dep.t.test.shrubland_desert,
                               dep.t.test.shrubland_temperate.forest,
                               dep.t.test.shrubland_boreal,
                               dep.t.test.desert_temperate.forest,
                               dep.t.test.desert_boreal,
                               dep.t.test.temperate_boreal)

pairwise.t.test.EDisp <- pairwise.t.test.EDisp[,1:3]

## ERich
## habitats
ed.t.test.close_open <- t.test(ed.close[,1], ed.open[,1])
## climates
ed.t.test.tropical_arid <- t.test(ed.tropical.climate[,1], ed.arid.climate[,1])
ed.t.test.tropical_temperate <- t.test(ed.tropical.climate[,1], ed.temperate.climate[,1])
ed.t.test.tropical_cold <- t.test(ed.tropical.climate[,1], ed.cold.climate[,1])
ed.t.test.arid_temperate <- t.test(ed.arid.climate[,1], ed.temperate.climate[,1])
ed.t.test.arid_cold <- t.test(ed.arid.climate[,1], ed.cold.climate[,1])
ed.t.test.temperate_cold <- t.test(ed.temperate.climate[,1], ed.cold.climate[,1])
### vegetations
ed.t.test.tropical.rainforest_tropical.seasonalforest <- t.test(ed.tropical.rainforest[,1], ed.tropical.seasonal.forest[,1])
ed.t.test.tropical.rainforest_savanna<-t.test(ed.tropical.rainforest[,1], ed.savanna[,1])
ed.t.test.tropical.rainforest_grassland <- t.test(ed.tropical.rainforest[,1], ed.grassland[,1])
ed.t.test.tropical.rainforest_shrubland <- t.test(ed.tropical.rainforest[,1], ed.shrubland[,1])
ed.t.test.tropical.rainforest_desert <- t.test(ed.tropical.rainforest[,1], ed.desert[,1])
ed.t.test.tropical.rainforest_temperate.forest <- t.test(ed.tropical.rainforest[,1], ed.temperate.forest[,1])
ed.t.test.tropical.rainforest_boreal <- t.test(ed.tropical.rainforest[,1], ed.boreal[,1])
ed.t.test.tropical.seasonalforest_savanna<- t.test(ed.tropical.seasonal.forest[,1], ed.savanna[,1])
ed.t.test.tropical.seasonalforest_grassland <- t.test(ed.tropical.seasonal.forest[,1], ed.grassland[,1])
ed.t.test.tropical.seasonalforest_shrubland <- t.test(ed.tropical.seasonal.forest[,1], ed.shrubland[,1])
ed.t.test.tropical.seasonalforest_desert <- t.test(ed.tropical.seasonal.forest[,1], ed.desert[,1])
ed.t.test.tropical.seasonalforest_temperate.forest <- t.test(ed.tropical.seasonal.forest[,1], ed.temperate.forest[,1])
ed.t.test.tropical.seasonalforest_boreal <- t.test(ed.tropical.seasonal.forest[,1], ed.boreal[,1])
ed.t.test.savanna_grassland <- t.test(ed.savanna[,1], ed.grassland[,1])
ed.t.test.savanna_shrubland <- t.test(ed.savanna[,1], ed.shrubland[,1])
ed.t.test.savanna_desert <- t.test(ed.savanna[,1], ed.desert[,1])
ed.t.test.savanna_temperate.forest <- t.test(ed.savanna[,1], ed.temperate.forest[,1])
ed.t.test.savanna_boreal <- t.test(ed.savanna[,1], ed.boreal[,1])
ed.t.test.grassland_shrubland <- t.test(ed.grassland[,1], ed.shrubland[,1])
ed.t.test.grassland_desert <- t.test(ed.grassland[,1], ed.desert[,1])
ed.t.test.grassland_temperate.forest <- t.test(ed.grassland[,1], ed.temperate.forest[,1])
ed.t.test.grassland_boreal <- t.test(ed.grassland[,1], ed.boreal[,1])
ed.t.test.shrubland_desert <- t.test(ed.shrubland[,1], ed.desert[,1])
ed.t.test.shrubland_temperate.forest <- t.test(ed.shrubland[,1], ed.temperate.forest[,1])
ed.t.test.shrubland_boreal <- t.test(ed.shrubland[,1], ed.boreal[,1])
ed.t.test.desert_temperate.forest <- t.test(ed.desert[,1], ed.temperate.forest[,1])
ed.t.test.desert_boreal <- t.test(ed.desert[,1], ed.boreal[,1])
ed.t.test.temperate_boreal <- t.test(ed.temperate.forest[,1], ed.boreal[,1])

pairwise.t.test.ERich <- rbind(ed.t.test.close_open,
                               ed.t.test.tropical_arid,
                               ed.t.test.tropical_temperate,
                               ed.t.test.tropical_cold,
                               ed.t.test.arid_temperate,
                               ed.t.test.arid_cold,
                               ed.t.test.temperate_cold,
                               ed.t.test.tropical.rainforest_tropical.seasonalforest,
                               ed.t.test.tropical.rainforest_savanna,
                               ed.t.test.tropical.rainforest_grassland,
                               ed.t.test.tropical.rainforest_shrubland,
                               ed.t.test.tropical.rainforest_desert,
                               ed.t.test.tropical.rainforest_temperate.forest,
                               ed.t.test.tropical.rainforest_boreal,
                               ed.t.test.tropical.seasonalforest_savanna,
                               ed.t.test.tropical.seasonalforest_grassland,
                               ed.t.test.tropical.seasonalforest_shrubland,
                               ed.t.test.tropical.seasonalforest_desert,
                               ed.t.test.tropical.seasonalforest_temperate.forest,
                               ed.t.test.tropical.seasonalforest_boreal,
                               ed.t.test.savanna_grassland,
                               ed.t.test.savanna_shrubland,
                               ed.t.test.savanna_desert,
                               ed.t.test.savanna_temperate.forest,
                               ed.t.test.savanna_boreal,
                               ed.t.test.grassland_shrubland,
                               ed.t.test.grassland_desert,
                               ed.t.test.grassland_temperate.forest,
                               ed.t.test.grassland_boreal,
                               ed.t.test.shrubland_desert,
                               ed.t.test.shrubland_temperate.forest,
                               ed.t.test.shrubland_boreal,
                               ed.t.test.desert_temperate.forest,
                               ed.t.test.desert_boreal,
                               ed.t.test.temperate_boreal)

pairwise.t.test.ERich <- pairwise.t.test.ERich[,1:3]

write.csv(pairwise.t.test.EDisp, file = "EDisp_pairwise_t_test_statistic_24July2-17.csv")
write.csv(pairwise.t.test.ERich, file = "ERich_pairwise_t_test_statistic_24July2-17.csv")


################################################
#### Average species in each enviornment type
################################################

mean.species.func <- function(input.data) {
  no.mean.species <- mean(table(input.data$Community.No))
  return(no.mean.species)
}

## based on two habitats
no.species.open <- mean.species.func(mcm.open)
no.species.close <- mean.species.func(mcm.close)
## based on four climates
no.species.tropical.climate <- mean.species.func(mcm.tropical.climate)
no.species.arid.climate <-mean.species.func(mcm.arid.climate)
no.species.temperate.climate <- mean.species.func(mcm.temperate.climate)
no.species.cold.climate <- mean.species.func(mcm.cold.climate)
## based on eight vegetation
no.species.tropical.rainforest <- mean.species.func(mcm.tropical.rainforest)
no.species.tropical.seasonal.forest <- mean.species.func(mcm.tropical.seasonal.forest)
no.species.savanna <- mean.species.func(mcm.savanna)
no.species.grassland <- mean.species.func(mcm.grassland)
no.species.shrubland <- mean.species.func(mcm.shrubland)
no.species.desert <- mean.species.func(mcm.desert)
no.species.temperate.forest <- mean.species.func(mcm.temperate.forest)
no.species.boreal.forest <- mean.species.func(mcm.boreal.forest)


##################################################
### Chinese Mesozoic mammal community analysis ###
##################################################

# Tiaojishan(TJS) community refers to Daxishan (DXS); Jiulongshan(JLS) community refers to Daohugou (DHG); based on
# Meng et al., 2015 (Arboreal docodont), adding Messel fauna
setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis")

ChMeMaCo<-read.csv(file="ChinaMesozoicMammalCommunity_28June2018.csv", header = T, sep=",")
ChMeMaCo$Eco.type <- with(ChMeMaCo, paste0(BSRank, DietRank, LocomotorRank))
length(unique(ChMeMaCo$Eco.type))

unique.ecotype.extinct <- unique(ChMeMaCo$Eco.type)
unique.ecotype.extinct <- data.frame(unique.ecotype.extinct)
colnames(unique.ecotype.extinct) <- c("ecotype")
unique.ecotype.extinct.splitted <- t(sapply(unique.ecotype.extinct$ecotype, 
                                            function(x) substring (x, 1:3, 1:3)))
unique.ecotype.extinct.splitted <- data.frame(unique.ecotype.extinct.splitted)
write.csv(unique.ecotype.extinct.splitted, file = "Unique_ecotype_extinct_20Sept2018.csv")


### different communities
data.dwzz<-ChMeMaCo[ChMeMaCo$Community=="DWZZ",][,4:6]   # DWZZ
write.table(data.dwzz, file="data.dwzz.csv", sep=",",row.names=F, col.names = F)
data.ljt_jsg <- ChMeMaCo[ChMeMaCo$Community=="LJT-JSG",][,4:6]   # LJT-JSG
write.table(data.ljt_jsg, file="data.ljt_jsg.csv", sep=",",row.names=F, col.names = F)
data.tjs<- ChMeMaCo[ChMeMaCo$Community=="TJS",][,4:6] # TJS
write.table(data.tjs, file="data.tjs.csv", sep=",", row.names=F, col.names = F)
data.jls<- ChMeMaCo[ChMeMaCo$Community=="JLS",][,4:6] #JLS
write.table(data.jls, file="data.jls.csv", sep=",", row.names=F, col.names = F)
data.msl<- ChMeMaCo[ChMeMaCo$Community=="MESL",][,4:6] # Messel
write.table(data.msl, file="data.msl.csv", sep=",", row.names=F, col.names = F)

# Ecological disparity function for extinct mammalian communities
data.eds<-data.frame()
Eds.extinct<- function (data.input) {
  for (i in 1:(dim(data.input)[1]-1)) { ### select first species to penultimate species within the community
    for (j in i:(dim(data.input)[1]-1)) {### select second to the last species with the communitiy
      data.eds.abs<-abs(data.input[i,]-data.input[j+1,]) ## absoluate different betweeen each function trait
      data.eds.pair<-rowSums(data.eds.abs) ### ecological disparity between a pair of species
      data.eds<-rbind(data.eds,c(data.eds.pair)) ### adding all results of all communities 
    }
  }
  return(data.eds)
}

# ecological disparity of each community
eds.dwzz<-Eds.extinct(data.dwzz)
eds.ljt_jsg<-Eds.extinct(data.ljt_jsg)
eds.tjs<-Eds.extinct(data.tjs)
eds.jls<-Eds.extinct(data.jls)
eds.msl<-Eds.extinct(data.msl)

## t.test among those commmunties
dwzz.ljt_jsg <- t.test(eds.dwzz, eds.ljt_jsg)
dwzz.tjs <- t.test(eds.dwzz, eds.tjs)
dwzz.jls <- t.test(eds.dwzz, eds.jls)
dwzz.msl <- t.test(eds.dwzz, eds.msl)
ljt_jsg.tjs <- t.test(eds.ljt_jsg, eds.tjs)
ljt_jsg.jls <- t.test(eds.ljt_jsg, eds.jls)
ljt_jsg.msl <- t.test(eds.ljt_jsg, eds.msl)
tjs.jls <- t.test(eds.tjs, eds.jls)
tjs.msl <- t.test(eds.tjs, eds.msl)
jls.msl <- t.test(eds.jls, eds.msl)

t.test.extinct <- rbind(dwzz.ljt_jsg,
                        dwzz.tjs,
                        dwzz.jls,
                        dwzz.msl,
                        ljt_jsg.tjs,
                        ljt_jsg.jls,
                        ljt_jsg.msl,
                        tjs.jls,
                        tjs.msl,
                        jls.msl)[,1:3]
colnames(t.test.extinct) <- c("t", "df", "P-value")
rownames(t.test.extinct) <- c("DWZZ vs LJT-JSG",
                              "DWZZ vs TJS",
                              "DWZZ vs JLS",
                              "DWZZ vs MSL",
                              "LJT-JSG vs TJS",
                              "LJT-JSG vs JLS",
                              "LJT-JSG vs MSL",
                              "TJS vs JLS",
                              "TJS vs MSL",
                              "JLS vs MSL")

write.csv(t.test.extinct, file = "t_test_extinct_communities_17Jan2018.csv")

# Ecological diversity function for extinct mammalian communities
data.edv<-data.frame()
Edv.extinct<-function (data.input) {
  data.edv.temp<-unique(data.input)
  data.edv<-rbind(data.edv,c(nrow(data.edv.temp)))
  return (data.edv)
}

# ecological diversity of each community
edv.dwzz<-Edv.extinct(data.dwzz)
edv.ljt_jsg<-Edv.extinct(data.ljt_jsg)
edv.tjs<-Edv.extinct(data.tjs)
edv.jls<-Edv.extinct(data.jls)
edv.msl<-Edv.extinct(data.msl)

###
# plots of ecological disparity and diversity of four extinct mammalian communities

par(mfcol=c(1,1),
    oma=c(1,1,0,1), 
    mar=c(1,2,0,1), 
    mgp=c(1,0.6,0),
    tck=-0.02)


mean.eds.extinct <- t(matrix(cbind(mean(eds.dwzz[,1]),
                                   mean(eds.ljt_jsg[,1]),
                                   mean(eds.tjs[,1]),
                                   mean(eds.jls[,1]),
                                   mean(eds.msl[,1])),nrow=5))
sd.eds.extinct <- t(matrix(cbind(sd(eds.dwzz[,1]),
                                 sd(eds.ljt_jsg[,1]),
                                 sd(eds.tjs[,1]),
                                 sd(eds.jls[,1]),
                                 sd(eds.msl[,1])),nrow=5))
# NO MEAN diversity for each extinct community; all together
all.edv.extinct<-t(matrix(cbind(edv.dwzz[,1],
                                edv.ljt_jsg[,1],
                                edv.tjs[,1],
                                edv.jls[,1],
                                edv.msl[,1])))
# length of ecological dispairity and diversity
length.eds.extinct<-t(matrix(cbind(length(eds.dwzz[,1]),
                                   length(eds.ljt_jsg[,1]),
                                   length(eds.tjs[,1]),
                                   length(eds.msl[,1])),nrow=5))

mean.eds.edv.extinct<-rbind(mean.eds.extinct, all.edv.extinct)
sd.eds.edv.extinct<-rbind(sd.eds.extinct,c(0,0,0,0,0))
length.eds.edv.extinct<-rbind(length.eds.extinct, c(1,1,1,1,1))
barx <- barplot(mean.eds.edv.extinct, beside=T,names.arg=c("DWZZ","LJT-JSG","TJS","JLS","MESL"),
                ylim=c(0,20), col=c("gold3","cadetblue"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,11))
## plot the error bar
error.bar(barx, mean.eds.edv.extinct, 1.96*sd.eds.edv.extinct/sqrt(length.eds.edv.extinct))

legend(4, 13.5,legend=c("Ecological disparity", "Ecological diversity"),
       fill=c("gold3","cadetblue"),horiz=F,cex=1.2,
       box.col = "white",border="white")
mtext(outer=F, side=3, line=-3, text="Extinct mammal communities", cex=1.5)


##################################################
#####  Single Edsip and ERich examples   #########
##################################################

### using the EDisp fucntion of extinct mammalian communities

### Kinabalu data
Kinabalu <- mcm.habitat[mcm.habitat$Community.No==46,]
write.csv(Kinabalu, file = "Data_Kinabalu.csv.csv")

Eco.type.data.Kinabalu <- as.data.frame(with(Kinabalu, paste0(BS.Rank, Diet.Rank, Locomotor.Rank)))
colnames(Eco.type.data.Kinabalu) <- c("EcoType")
Eco.type.data.Kinabalu <-cbind(Kinabalu$Community.No, Kinabalu$Habitat, Kinabalu$Climate,
                               Kinabalu$Vegetation, Eco.type.data.Kinabalu)
colnames(Eco.type.data.Kinabalu) <- c("Community", "Habitat", "Climate", "Vegetation", "Ecotype")

### EDisp
Kinabalu.EDisp<-Eds.extinct(Kinabalu[,c(3:5)])
Kinabalu.EDisp.mean <- mean(Kinabalu.EDisp[,1])
Kinabalu.EDisp.sd <- sd(Kinabalu.EDisp[,1])
### ERich
Kinabalu.ERich <- length(unique(Eco.type.data.Kinabalu$Ecotype))


### Chihuahuan data
Chihuahuan<-mcm.habitat[mcm.habitat$Community.No==93,]
write.csv(Chihuahuan, file = "Data_Chihuahuan.csv.csv")

Eco.type.data.Chihuahuan <- as.data.frame(with(Chihuahuan, paste0(BS.Rank, Diet.Rank, Locomotor.Rank)))
colnames(Eco.type.data.Chihuahuan) <- c("EcoType")
Eco.type.data.Chihuahuan <-cbind(Chihuahuan$Community.No, Chihuahuan$Habitat, Chihuahuan$Climate,
                                 Chihuahuan$Vegetation, Eco.type.data.Chihuahuan)
colnames(Eco.type.data.Chihuahuan) <- c("Community", "Habitat", "Climate", "Vegetation", "Ecotype")

### EDisp
Chihuahuan.EDisp<-Eds.extinct(Chihuahuan[,c(3:5)])
Chihuahuan.EDisp.mean <- mean(Chihuahuan.EDisp[,1])
Chihuahuan.EDisp.sd <- sd(Chihuahuan.EDisp[,1])
### ERich
Chihuahuan.ERich <- length(unique(Eco.type.data.Chihuahuan$Ecotype))



