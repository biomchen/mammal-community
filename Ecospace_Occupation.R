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

setwd("/Users/mengchen/Documents/Research Projects/MM Community Project/Database and Analysis")

###########################################
# investigate the bs relationships with environmental factors
###################################################
bs_enviro_data<-read.csv(file="BS_vs_Enviro_factors_19022016.csv", sep = ",", header = T)
summary(lm(BS~abs(Lat), data = bs_enviro_data))
summary(lm(BS~Lon, data = bs_enviro_data)) # this linear relationship is becaause the data sample along longitude
summary(lm(BS~MAP, data = bs_enviro_data))
summary(lm(BS~MAT, data = bs_enviro_data))
summary(lm(BS~Elevation, data = bs_enviro_data))

### Read the database of the 45 communities
mcm<-read.csv(file="EcoDataOfAllComm_11Jan2016.csv", sep=",", header=T)

## reorder data based on the latitudinal gradients
mcm.reordered<-mcm[order(mcm$ClimateCode),]

# ### bootstrap about ecological disparity and diversity using funcitons
# source("EcoStructFunc_Bootstrap.R")
# results<-EcoFunc.eds.bs(mcm.reordered,340,1000)


## Write the recordered data
write.table(mcm.reordered,file="EcoDataOfAll_reordered.csv", sep=",",row.names=F,col.names=F)

mcm.data.reordered<-mcm.reordered ## because previous usage, so does here
cc<-unique(mcm.reordered$ClimateCode) ## climate code

##### The data of each climate region
mcm.reordered.tropical<-mcm.reordered[mcm.reordered$ClimateCode < 200,]
mcm.reordered.arid<-mcm.reordered[(mcm.reordered$ClimateCode >200) & 
																 	(mcm.reordered$ClimateCode < 300),]
mcm.reordered.temperate<-mcm.reordered[(mcm.reordered$ClimateCode >300) & 
																				(mcm.reordered$ClimateCode < 400),]
mcm.reordered.cold<-mcm.reordered[mcm.reordered$ClimateCode >400,]

# summary each rank of each function has how many occurrences
# body size rank
table(mcm.reordered.tropical$BodySizeRank)
table(mcm.reordered.arid$BodySizeRank)
table(mcm.reordered.temperate$BodySizeRank)
table(mcm.reordered.cold$BodySizeRank)
# diet rank
table(mcm.reordered.tropical$DietRank)
table(mcm.reordered.arid$DietRank)
table(mcm.reordered.temperate$DietRank)
table(mcm.reordered.cold$DietRank)
# locomotor rank
table(mcm.reordered.tropical$LocomotorRank)
table(mcm.reordered.arid$LocomotorRank)
table(mcm.reordered.temperate$LocomotorRank)
table(mcm.reordered.cold$LocomotorRank)

# Chi-Square test for
chisq.table.bodysize<-matrix(c(14,76,96,25,36,  ### Tropical  Those data from table(mcm.reordered.tropical[,3])
															 24,20,11,7,8,      ### Arid
															 75,61,26,14,9,   ### Temperate
															 13,17,9,15,13),  ### Cold
														 nrow=4,ncol=5)

chisq.table.diet<-matrix(c(17,42,117,14,28,29,  ### Tropical Those data from table(mcm.reordered.tropical[,4])
													  17,15,22,8,0,8,       ### Arid
														16,26,72,18,2,51,   ### Temperate
														14,9,22,1,1,20),  ### Cold
														 nrow=4,ncol=6)
chisq.table.locomotion<-matrix(c(2,51,60,124,4,3,2,1,  ### Tropical Those data from table(mcm.reordered.tropical[,5])
																 0,0,14,34,0,11,3,8,       ### Arid
																 1,5,45,85,5,20,12,12,   ### Temperate
																 4,0,16,41,5,0,1,0),  ### Cold
															 nrow=4,ncol=8)

result.chisq.bodysize<-chisq.test(chisq.table.bodysize)
result.chisq.diet<-chisq.test(chisq.table.diet)
result.chisq.locomotion<-chisq.test(chisq.table.locomotion)

#################################################
## Ecological Disparity and Diversity Analysis ##
#################################################


## Global disparity and diversity
ed.global<-data.frame()
dep.global<-data.frame()
cn<-unique(mcm.data.reordered$Community) ## The number of commuities
ed.global.data<-data.frame()
mean.dep.community<-data.frame()
each.community.mean.dep.global<-data.frame()


for (j in cn) {
	newdata.global<-mcm.data.reordered[mcm.data.reordered$Community==j,]
	ed.global.data<-rbind(ed.global.data,c(newdata.global[c(3,4,5)])) ### entire ecological raw dataset
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
	    dep.global<-rbind(dep.global,c(dep.temp.global))
	  }
	}
}



#### global Mean ecological disparity
mean.dep.global <- mean (dep.global[,1])

### the true ecological occupations
unique.global.data<-ed.global[!duplicated(ed.global),]
### mean value of duplicates
mean.global.duplicates<-nrow(ed.global)/nrow(unique.global.data)

###################################################################################################
##### relationship test between three ecological paramters using kendall's concordance analysis
kendall.global(ed.global.data, nperm = 999, mult = "holm")
kendall.post(ed.global.data, nperm = 999, mult = "holm")
cor(ed.global.data, method="kendall", use="pairwise")

### Tropical; climate code start with 1XX
ed.tropical<-data.frame()
dep.tropical<-data.frame()
mcm.tropical<-data.frame()
mcm.all.tropical<-data.frame()
ed.tropical.data<-data.frame()


for (e in 1:3) {
	mcm.temp.tropical<-mcm.reordered*(mcm.reordered$ClimateCode==cc[e])
	mcm.sub.tropical<-subset (mcm.temp.tropical, ClimateCode > 0) ## Remove those 0,0,0 rows
	mcm.all.tropical<-rbind(mcm.all.tropical,c(mcm.sub.tropical))
	n.tropical<-unique(mcm.all.tropical$Community)
}
	
for (j in 1:length(n.tropical)) {
	mcm.all.temp.tropical<-mcm.all.tropical*(mcm.all.tropical$Community==n.tropical[j])
	mcm.community.temp.tropical<-subset(mcm.all.temp.tropical,Community > 0)
	### ecological raw dataset
	ed.tropical.data<-rbind(ed.tropical.data,c(mcm.community.temp.tropical[,c(3,4,5)]))
	## diversity
	data.tropical<-mcm.community.temp.tropical[,c(3,4,5)]
	result.tropical<-nrow(unique(data.tropical))
	ed.tropical<-rbind(ed.tropical,c(result.tropical)) 
	#### disparity
	for (u in 1:(dim(data.tropical)[1]-1)) {
	  for (v in u:(dim(data.tropical)[1]-1)) {
	    results.temp.tropical<-abs(data.tropical[u,]-data.tropical[v+1,])
	    dep.temp.tropical<-rowSums(results.temp.tropical)
	    dep.tropical<-rbind(dep.tropical,c(dep.temp.tropical))
	  }
	}
}


### the true ecological occupations
unique.tropical.data<-ed.tropical.data[!duplicated(ed.tropical.data),]
### mean value of duplicates
mean.tropical.duplicates<-nrow(ed.tropical.data)/nrow(unique.tropical.data)

### Arid; climate code start with 2XX
ed.arid<-data.frame()
dep.arid<-data.frame()
mcm.arid<-data.frame()
mcm.all.arid<-data.frame()
ed.arid.data<-data.frame()

for (e in 4:7) {
	mcm.temp.arid<-mcm.reordered*(mcm.reordered$ClimateCode==cc[e])
	mcm.sub.arid<-subset (mcm.temp.arid, ClimateCode > 0) ## Remove those 0,0,0 rows
	mcm.all.arid<-rbind(mcm.all.arid,c(mcm.sub.arid))
	n.arid<-unique(mcm.all.arid$Community)
}	

for (g in 1:length(n.arid)) {
	mcm.all.temp.arid<-mcm.all.arid*(mcm.all.arid$Community==n.arid[g])
	mcm.community.temp.arid<-subset(mcm.all.temp.arid,Community > 0)
	### ecological raw dataset
	ed.arid.data<-rbind(ed.arid.data,c(mcm.community.temp.arid[,c(3,4,5)]))
	## diversity
	data.arid<-mcm.community.temp.arid[,c(3,4,5)]
	result.arid<-nrow(unique(data.arid))
	ed.arid<-rbind(ed.arid,c(result.arid)) 
	#### disparity
	for (u in 1:(dim(data.arid)[1]-1)) {
	  for (v in u:(dim(data.arid)[1]-1)) {
	    results.temp.arid<-abs(data.arid[u,]-data.arid[v+1,])
			dep.temp.arid<-rowSums(results.temp.arid)
			dep.arid<-rbind(dep.arid,c(dep.temp.arid))
		}
	}
}
### the true ecological occupations
unique.arid.data<-ed.arid.data[!duplicated(ed.arid.data),]
### mean value of duplicates
mean.arid.duplicates<-nrow(ed.arid.data)/nrow(unique.arid.data)

		 
### Temperate; climate code start with 3XX
ed.temperate<-data.frame()
dep.temperate<-data.frame()
mcm.temperate<-data.frame()
mcm.all.temperate<-data.frame()
ed.temperate.data<-data.frame()

for (e in 8:14) {
	mcm.temp.temperate<-mcm.reordered*(mcm.reordered$ClimateCode==cc[e])
	mcm.sub.temperate<-subset (mcm.temp.temperate, ClimateCode > 0) ## Remove those 0,0,0 rows
	mcm.all.temperate<-rbind(mcm.all.temperate,c(mcm.sub.temperate))
	n.temperate<-unique(mcm.all.temperate$Community)
}

for (h in 1:length(n.temperate)) {
	mcm.all.temp.temperate<-mcm.all.temperate*(mcm.all.temperate$Community==n.temperate[h])
	mcm.community.temp.temperate<-subset(mcm.all.temp.temperate,Community > 0)
	### diversity raw data
	ed.temperate.data<-rbind(ed.temperate.data,c(mcm.community.temp.temperate[,c(3,4,5)]))
	data.temperate<-mcm.community.temp.temperate[,c(3,4,5)]
	result.temperate<-nrow(unique(data.temperate))
	## diversity
	ed.temperate<-rbind(ed.temperate,c(result.temperate)) 
	#### disparity
	for (u in 1:(dim(data.temperate)[1]-1)) {
		for (v in u:(dim(data.temperate)[1]-1)) {
			results.temp.temperate<-abs(data.temperate[u,]-data.temperate[v+1,])
			dep.temp.temperate<-rowSums(results.temp.temperate)
			dep.temperate<-rbind(dep.temperate,c(dep.temp.temperate))
		}
	}
}

### the true ecological occupations
unique.temperate.data<-ed.temperate.data[!duplicated(ed.temperate.data),]
### mean value of duplicates
mean.temperate.duplicates<-nrow(ed.temperate.data)/nrow(unique.temperate.data)


### Cold; climate code start with 4XX
ed.cold<-data.frame()
dep.cold<-data.frame()
mcm.cold<-data.frame()
mcm.all.cold<-data.frame()
ed.cold.data<-data.frame()

for (e in 15:17) {
	mcm.temp.cold<-mcm.reordered*(mcm.reordered$ClimateCode==cc[e])
	mcm.sub.cold<-subset (mcm.temp.cold, ClimateCode > 0)## Remove those 0,0,0 rows
	mcm.all.cold<-rbind(mcm.all.cold,c(mcm.sub.cold))
	n.cold<-unique(mcm.all.cold$Community)
}

for (h in 1:length(n.cold)) {
	mcm.all.temp.cold<-mcm.all.cold*(mcm.all.cold$Community==n.cold[h])
	mcm.community.temp.cold<-subset(mcm.all.temp.cold,Community > 0)
	### ecological raw dataset
	ed.cold.data<-rbind(ed.cold.data,c(mcm.community.temp.cold[,c(3,4,5)]))
	data.cold<-mcm.community.temp.cold[,c(3,4,5)]
	result.cold<-nrow(unique(data.cold))
	## diversity
	ed.cold<-rbind(ed.cold,c(result.cold))
	#### disparity
	for (u in 1:(dim(data.cold)[1]-1)) {
		for (v in u:(dim(data.cold)[1]-1)) {
			results.temp.cold<-abs(data.cold[u,]-data.cold[v+1,])
			dep.temp.cold<-rowSums(results.temp.cold)
			dep.cold<-rbind(dep.cold,c(dep.temp.cold))
		}
	}
}

### the true ecological occupations
unique.cold.data<-ed.cold.data[!duplicated(ed.cold.data),]
### mean value of duplicates
mean.cold.duplicates<-nrow(ed.cold.data)/nrow(unique.cold.data)

##################################################
# plot the ecological parameters using percentage
#################################################
setwd("/Users/mengchen/Documents/Research Projects/MM Community Project/Database and Analysis")

eco.para.data<-read.csv(file = "EcoParameter_Percent_14022016.csv", sep = "," , header = T)

par(mfcol=c(1,1),
        oma=c(0,0,0,0), 
        mar=c(3,3,1,1), 
        mgp=c(1,0.6,0),
        tck=-0.02)
# Circles is vectors all radii of circiles, the radii should be put together as 1 in total
## Tropical  
symbols(x=1:8, y=c(1,1,1,1,1,1,1,1), circles=sqrt(eco.para.data[1,13:20]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="orange", inches = F)
symbols(x=1:6, y=c(2,2,2,2,2,2), circles=sqrt(eco.para.data[1,7:12]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="orange", add = T, inches = F)
symbols(x=1:5, y=c(3,3,3,3,3), circles=sqrt(eco.para.data[1,2:6]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="orange", add = T, inches = F)
## arid
symbols(x=1:8, y=c(1,1,1,1,1,1,1,1), circles=sqrt(eco.para.data[2,13:20]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="grey60", inches = F)
symbols(x=1:6, y=c(2,2,2,2,2,2), circles=sqrt(eco.para.data[2,7:12]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="grey60", add = T, inches = F)
symbols(x=1:5, y=c(3,3,3,3,3), circles=sqrt(eco.para.data[2,2:6]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="grey60", add = T, inches = F)
# temperate
symbols(x=1:8, y=c(1,1,1,1,1,1,1,1), circles=sqrt(eco.para.data[3,13:20]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="green3", inches = F)
symbols(x=1:6, y=c(2,2,2,2,2,2), circles=sqrt(eco.para.data[3,7:12]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="green3", add = T, inches = F)
symbols(x=1:5, y=c(3,3,3,3,3), circles=sqrt(eco.para.data[3,2:6]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="green3", add = T, inches = F)
# cold
symbols(x=1:8, y=c(1,1,1,1,1,1,1,1), circles=sqrt(eco.para.data[4,13:20]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="dodgerblue3", inches = F)
symbols(x=1:6, y=c(2,2,2,2,2,2), circles=sqrt(eco.para.data[4,7:12]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="dodgerblue3", add = T, inches = F)
symbols(x=1:5, y=c(3,3,3,3,3), circles=sqrt(eco.para.data[4,2:6]/pi), xlim=c(0,8), ylim=c(0,4), ann=F, bg="dodgerblue3", add = T, inches = F)
# legend; 0.1+0.2+0.3+0.4=1
symbols(x=1:6, y=c(1,1,1,1,1,1), circles=sqrt(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)/pi), xlim=c(0,8), ylim=c(0,2), ann=F, bg="white", inches = F)

###########################################################################################
# investigate what kind of ecospace has been occupied by SMC in four different environments
###########################################################################################
# tropical
eco.func.tropical<-read.csv(file="EcoFuncType_tropical_15022016.csv", sep = "," , header = F)
length(table(eco.func.tropical))
# arid
eco.func.arid<-read.csv(file="EcoFuncType_arid_15022016.csv", sep = "," , header = F)
length(table(eco.func.arid))
# temperate
eco.func.temperate<-read.csv(file="EcoFuncType_temperate_15022016.csv", sep = "," , header = F)
length(table(eco.func.temperate))
# cold
eco.func.cold<-read.csv(file="EcoFuncType_cold_15022016.csv", sep = "," , header = F)
length(table(eco.func.cold))


#### Save the body size, diet, and locomotion occurances in the dataset

write.table(ed.global.data, file="EdGlobalData.txt",sep="")
write.table(ed.tropical.data, file = "EdTropicalData.txt",sep="")
write.table(ed.tropical.data,file="EdTropicalDataAll.csv", sep=",",row.names=F,col.names=F)
write.table(ed.arid.data, file = "EdAridData.txt",sep="")
write.table(ed.arid.data,file="EdAridDataAll.csv", sep=",",row.names=F,col.names=F)
write.table(ed.temperate.data, file = "EdTemperateData.txt",sep="")
write.table(ed.temperate.data,file="EdTemperateDataAll.csv", sep=",",row.names=F,col.names=F)
write.table(ed.cold.data, file = "EdColdData.txt",sep="")
write.table(ed.cold.data,file="EdColdDataAll.csv", sep=",",row.names=F,col.names=F)

### The saved data files are transfered into the single ecoloical combinations
### manually using Excel import functions

#### Density plot of the ecological diversity of each region

par(mfcol=c(4,1),
		oma=c(3,3,1,2), 
		mar=c(2,2,1,1), 
		mgp=c(1,0.6,0),
		tck=-0.02)

## Tropical
plot(density(ed.tropical.data[,3],bw=0.4),xlim=c(1,8),
		 ylim=c(0,0.6),col="skyblue",las=1,lwd=1,bty="n",
		 main="",xlab="",ylab="",xaxt="n")
lines(density(ed.tropical.data[,2],bw=0.4),col="orange")
lines(density(ed.tropical.data[,1],bw=0.4),col="red")
mtext(outer=F,side=3,line=-1,text="Tropical",at=7,cex=0.8)
## Arid
plot(density(ed.arid.data[,3],bw=0.4),xlim=c(1,8),
		 ylim=c(0,0.6),col="skyblue",las=1,lwd=1,bty="n",
		 main="",xlab="",ylab="",xaxt="n")
lines(density(ed.arid.data[,2],bw=0.4),col="orange")
lines(density(ed.arid.data[,1],bw=0.4),col="red")
mtext(outer=F,side=3,line=-1,text="Arid",at=7,cex=0.8)
## Temperate
plot(density(ed.temperate.data[,3],bw=0.4),xlim=c(1,8),
		 ylim=c(0,0.6),col="skyblue",las=1,lwd=1,bty="n",
		 main="",xlab="",ylab="",xaxt="n")
lines(density(ed.temperate.data[,2],bw=0.4),col="orange")
lines(density(ed.temperate.data[,1],bw=0.4),col="red")
mtext(outer=F,side=3,line=-1,text="Temperate",at=7,cex=0.8)
## Cold
plot(density(ed.cold.data[,3],bw=0.4),xlim=c(1,8),
		 ylim=c(0,0.6),col="skyblue",las=1,lwd=1,bty="n",
		 main="",xlab="",ylab="")
lines(density(ed.cold.data[,2],bw=0.4),col="orange")
lines(density(ed.cold.data[,1],bw=0.4),col="red")
mtext(outer=F,side=3,line=-1,text="Cold",at=7,cex=0.8)

mtext(outer=T,side=1,line=0.5,text="Ecological Parameter Rank",cex=0.8)
mtext(outer=T,side=2,line=1,text="Density",at=0.5,cex=0.8)

#####################################################
# Plot of diversity based on different climate zone
#####################################################
par(mfcol=c(1,1),
    oma=c(3,3,1,2), 
    mar=c(4,2,1,1), 
    mgp=c(1,0.6,0),
    tck=-0.02)

# color coding for four environments
col.global1<-c(rep("#f1595f", 15),  # number indicate how many communities in that environment
               rep("#79c36a", 6),
               rep("#599ad3", 13),
               rep("#f9a65a", 6))

community.no<-unique(mcm.reordered$Community)

barplot(ed.global[,1],las=1, names.arg=community.no,col=col.global1)

mtext(outer=F,side=1,line=2.5,text="Tropical", col="#f1595f", at=9.5, font=2, cex=1.5)
mtext(outer=F,side=1,line=2.5,text="Arid", col="#79c36a", at=21.5, font=2, cex=1.5)
mtext(outer=F,side=1,line=2.5,text="Temperate", col="#599ad3", at=34, font=2, cex=1.5)
mtext(outer=F,side=1,line=2.5,text="Cold", col="#f9a65a", at=43.8, font=2, cex=1.5)

mtext(outer=T,side=2,line=0.5,text="Ecological Diversity",las=3, font=2, cex=1.5)

# MEAN DISPARITY
# Global mean disparity
mean.global<-data.frame()
sd.global<-data.frame()

mean.temp.global<-apply(as.matrix(dep.global[,1]),2,mean)
mean.global<-rbind(mean.global,c(mean.temp.global))
sd.temp.global<-apply(as.matrix(dep.global[,1]),2,sd)
sd.global<-rbind(sd.global,c(sd.temp.global))

mean.global<-t(mean.global)
sd.global<-t(sd.global)

## Tropical Mean Disparity
mean.tropical<-data.frame()
sd.tropical<-data.frame()

mean.temp.tropical<-apply(as.matrix(dep.tropical[,1]),2,mean)
mean.tropical<-rbind(mean.tropical,c(mean.temp.tropical))
sd.temp.tropical<-apply(as.matrix(dep.tropical[,1]),2,sd)
sd.tropical<-rbind(sd.tropical,c(sd.temp.tropical))

mean.tropical<-t(mean.tropical)
sd.tropical<-t(sd.tropical)

## Arid Mean Disparity
mean.arid<-data.frame()
sd.arid<-data.frame()

mean.temp.arid<-apply(as.matrix(dep.arid[,1]),2,mean)
mean.arid<-rbind(mean.arid,c(mean.temp.arid))
sd.temp.arid<-apply(as.matrix(dep.arid[,1]),2,sd)
sd.arid<-rbind(sd.arid,c(sd.temp.arid))

mean.arid<-t(mean.arid)
sd.arid<-t(sd.arid)

## Temperate Mean Disparity
mean.temperate<-data.frame()
sd.temperate<-data.frame()

mean.temp.temperate<-apply(as.matrix(dep.temperate[,1]),2,mean)
mean.temperate<-rbind(mean.temperate,c(mean.temp.temperate))
sd.temp.temperate<-apply(as.matrix(dep.temperate[,1]),2,sd)
sd.temperate<-rbind(sd.temperate,c(sd.temp.temperate))

mean.temperate<-t(mean.temperate)
sd.temperate<-t(sd.temperate)

## Cold Mean Disparity
mean.cold<-data.frame()
sd.cold<-data.frame()

mean.temp.cold<-apply(as.matrix(dep.cold[,1]),2,mean)
mean.cold<-rbind(mean.cold,c(mean.temp.cold))
sd.temp.cold<-apply(as.matrix(dep.cold[,1]),2,sd)
sd.cold<-rbind(sd.cold,c(sd.temp.cold))

mean.cold<-t(mean.cold)
sd.cold<-t(sd.cold)

## Barplot of Disparity
## funciton of standard error in barplot
error.bar <- function(x, y, upper, lower=upper, length=0.1,...) {
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

mean.all <- cbind(mean.global,mean.tropical,mean.arid,mean.temperate,mean.cold)
sd.all <- cbind(sd.global,sd.tropical,sd.arid,sd.temperate,sd.cold)
mean.names <- c("Global","Tropical","Arid","Temperate","Cold")
length.all <- cbind(length(t(dep.global)),length(t(dep.tropical)),
                       length(t(dep.arid)),length(t(dep.temperate)),length(t(dep.cold)))

barx <- barplot(mean.all, names.arg=mean.names,ylim=c(0,6), col="grey70", 
                cex.axis=0.8,axis.lty=1, space=1, cex.names=0.8,las=1)
## plot the error bar
error.bar(barx, mean.all, 1.96*sd.all/sqrt(length.all))

mtext(outer=F, side=2, line=1.5, text="Mean Disparity",cex=0.8)
mtext(outer=F, side=3, line=-2, at=5,cex=0.8)


#### Distribution of the disparity
par(mfcol=c(4,1),
		oma=c(3,3,0,0), 
		mar=c(2,2,2,1), 
		mgp=c(1,0.6,0),
		tck=-0.02,
		bty="n",
		yaxs="i")

##Tropical
h.tropical<-hist(dep.tropical[,1],xlim=c(0,14),ylim=c(0,600),breaks=14,
								 main="",xlab="",ylab="",las=1,col="grey50")
xfit.tropical<-seq(min(dep.tropical[,1]),max(dep.tropical[,1]),length=40) 
yfit.tropical<-dnorm(xfit.tropical,mean=mean(dep.tropical[,1]),sd=sd(dep.tropical[,1])) 
yfit.tropical <- yfit.tropical*diff(h.tropical$mids[1:2])*length(dep.tropical[,1]) 
lines(xfit.tropical, yfit.tropical, col="blue", lwd=1.5)
abline(v=median(dep.global[,1]),lty=5,lwd=1,col="grey50")
mtext(outer=F, side=3, line=-2, at=10,text="Tropical \nN=2436",cex=0.9)

##Arid
h.arid<-hist(dep.arid[,1],xlim=c(0,14),ylim=c(0,30),breaks=14,
						 main="",xlab="",ylab="",las=1,col="grey50")
xfit.arid<-seq(min(dep.arid[,1]),max(dep.arid[,1]),length=40) 
yfit.arid<-dnorm(xfit.arid,mean=mean(dep.arid[,1]),sd=sd(dep.arid[,1])) 
yfit.arid <- yfit.arid*diff(h.arid$mids[1:2])*length(dep.arid[1]) 
lines(xfit.arid, yfit.arid, col="blue", lwd=1.5)
abline(v=median(dep.global[,1]),lty=5,lwd=1,col="grey50")
mtext(outer=F, side=3, line=-2, at=10,text="Arid \nN=272",cex=0.9)

## Temperate
h.temperate<-hist(dep.temperate[,1],xlim=c(0,14),ylim=c(0,250),breaks=14,
									main="",xlab="",ylab="",las=1,col="grey50")
xfit.temperate<-seq(min(dep.temperate[,1]),max(dep.temperate[,1]),length=40) 
yfit.temperate<-dnorm(xfit.temperate,mean=mean(dep.temperate[,1]),sd=sd(dep.temperate[,1])) 
yfit.temperate <- yfit.temperate*diff(h.temperate$mids[1:2])*length(dep.temperate[,1]) 
lines(xfit.temperate, yfit.temperate, col="blue", lwd=1.5)
abline(v=median(dep.global[,1]),lty=5,lwd=1,col="grey50")
mtext(outer=F, side=3, line=-2, at=10,text="Temperate \nN=1277",cex=0.9)

## Cold
h.cold<-hist(dep.cold[,1],xlim=c(0,14),ylim=c(0,100),breaks=14,
						 main="",xlab="",ylab="",las=1,col="grey50")
xfit.cold<-seq(min(dep.cold[,1]),max(dep.cold[,1]),length=40) 
yfit.cold<-dnorm(xfit.cold,mean=mean(dep.cold[,1]),sd=sd(dep.cold[,1])) 
yfit.cold<- yfit.cold*diff(h.cold$mids[1:2])*length(dep.cold[,1]) 
lines(xfit.cold, yfit.cold, col="blue", lwd=1.5)
abline(v=median(dep.global[,1]),lty=5,lwd=1,col="grey50")
mtext(outer=F, side=3, line=-2, at=10,text="Cold \nN=477",cex=0.9)

mtext(outer=T, side=2, line=1, text="Frequency",cex=0.9)
mtext(outer=T, side=1, line=0.5, text="Ecological Disparity",cex=0.9)


##### Density plots ########
par(mfcol=c(1,1),
		oma=c(2,2,0,0), 
		mar=c(2,2,2,1), 
		mgp=c(1,0.6,0),
		tck=-0.015,
		bty="n",
		yaxs="i")

# Plot 1
plot(density(dep.tropical[,1],bw=0.7),main="",xlab="",ylab="", ylim=c(0,0.22),
		 xlim=c(0,max(dep.global[,1])),las=1,col="#f1595f",lwd=2, cex.axis=0.7)
lines(density(dep.arid[,1],bw=0.7),main="",ylab="",
			las=1,col="#79c36a",lwd=2)
lines(density(dep.temperate[,1],bw=0.7),main="",ylab="",
			las=1,col="#599ad3",lwd=2)
lines(density(dep.cold[,1],bw=0.7),main="",ylab="",
			las=1,col="#f9a65a",lwd=2)

abline(v=mean(dep.global[,1]),lty=2)

legend(8.5,0.2,c("Tropical, N= 1762",
								 "Arid, N= 127",
								 "Temperate, N= 1203",
								 "Cold, N= 449"),
			 lty=1, lwd=2, bty="n",
			 col=c("#f1595f","#79c36a","#599ad3","#f9a65a","black"),
			 cex=0.8)

mtext(outer=T, side=2, line=0.5, text="Estimate Probability",cex=1)
mtext(outer=T, side=1, line=0.5, text="Ecological Disparity",cex=1)


######################################################################
### Test if the ecological disparity is normal distribution (skewness)
### Similarly the ecological diversity can be test in the same way
######################################################################
library(e1071)
skewness<-c(skewness(dep.global[,1]),skewness(dep.tropical[,1]),
						skewness(dep.arid[,1]),skewness(dep.temperate[,1]),
						skewness(dep.cold[,1]))

kurtosis<-c(kurtosis(dep.global[,1]),kurtosis(dep.tropical[,1]),
						kurtosis(dep.arid[,1]),kurtosis(dep.temperate[,1]),
						kurtosis(dep.cold[,1]))

s.error<-c(sd.global[,1]/sqrt(length(dep.global[,1])),
					 sd.tropical[,1]/sqrt(length(dep.tropical[,1])),
					 sd.arid[,1]/sqrt(length(dep.arid[,1])),
					 sd.temperate[,1]/sqrt(length(dep.temperate[,1])),
					 sd.cold[,1]/sqrt(length(dep.cold[,1])))

s.deviation<-c(sd.global[,1],sd.tropical[,1],
							 sd.arid[,1],sd.temperate[,1],
							 sd.cold[,1])

### Results show the ecological disparity (as well diversity) is not normally distributed.


#####################
## MEAN Diversity ###
#####################

## Tropical Mean Disparity

mean.tropical.ed<-mean(ed.tropical[,1])
sd.tropical.ed<-sd(ed.tropical[,1])
mean.tropical.ed<-t(mean.tropical.ed)
sd.tropical.ed<-t(sd.tropical.ed)

## Arid Mean Diversity

mean.arid.ed<-mean(ed.arid[,1])
sd.arid.ed<-sd(ed.arid[,1])
mean.arid.ed<-t(mean.arid.ed)
sd.arid.ed<-t(sd.arid.ed)

## Temperate Mean Diversity

mean.temperate.ed<-mean(ed.temperate[,1])
sd.temperate.ed<-sd(ed.temperate[,1])
mean.temperate.ed<-t(mean.temperate.ed)
sd.temperate.ed<-t(sd.temperate.ed)

## Cold Mean Diversity

mean.cold.ed<-mean(ed.cold[,1])
sd.cold.ed<-sd(ed.cold[,1])
mean.cold.ed<-t(mean.cold.ed)
sd.cold.ed<-t(sd.cold.ed)


###############################
### plot Diversity together ###
###############################
mean.ed.all<-cbind(mean.tropical.ed[,1],
							 mean.arid.ed[,1],
							 mean.temperate.ed[,1],
							 mean.cold.ed[,1])

sd.ed.all<-cbind(sd.tropical.ed[,1],
						 sd.arid.ed[,1],
						 sd.temperate.ed[,1],
						 sd.cold.ed[,1])

length.ed.all<-t(matrix(cbind(length(ed.tropical[,1]),
                              length(ed.arid[,1]),length(ed.temperate[,1]),
                              length(ed.cold[,1]))))

## all together
barx <- barplot(mean.ed.all, names.arg=c("Tropical",
																			"Arid","Temperate","Cold"), 
								ylim=c(0,20), col="grey70", cex.axis=0.8,
								axis.lty=1, cex.names=0.8,
								las=1,beside=TRUE)
## plot the error bar
error.bar(barx, mean.ed.all, 1.96*sd.ed.all/sqrt(length.ed.all))
mtext(outer=F, side=2, line=1.5, text="Mean Diversity",cex=0.8)


#################################################################################################
# The CA analysis for niche theory (Braak 1985): functional trait combination along environmental
# gradients (Latitudinal, Longtitudinal, elevational, MAT, MAP or etc)
#################################################################################################
community.ed.global.data<-cbind(mcm.data.reordered$Community,ed.global.data)
# write.table(community.ed.global.data, file="GlobalCommunityEdData_03Feb2016.csv", sep=",")

# there is a transformation of the data using EXCEL to combination functional trait cells
Community.FuncType<-read.csv(file="GlobalCommunityEdData_03Feb2016.csv", sep=",", header = T)


ClimateCode<-mcm.data.reordered$ClimateCode
Community.FuncType.ClimateCode<-cbind(Community.FuncType,
                                      ClimateCode)
Community.FuncType.Matrix<-table(Community.FuncType)

# then reorganize the Community.FuncType.Matrix
Community.FuncType.Matrix.reordered<-Community.FuncType.Matrix[unique(mcm.data.reordered$Community),]
Community.order<-unique(Community.FuncType$No)

jaccard.dis.global<-data.frame()
first.func.community<-data.frame()
second.func.community<-data.frame()

for (g in 1:(length(Community.order)-1)) {
  for (h in g:(length(Community.order)-1)) {
    first.func.community<-unique(Community.FuncType[Community.FuncType$No==Community.order[g],]$FuncType)
    second.func.community<-unique(Community.FuncType[Community.FuncType$No==Community.order[h+1],]$FuncType)
    # jaccarad disimiliarity equation
    jaccard.dis.temp<-1-length(intersect(first.func.community,second.func.community))/length(union(first.func.community,second.func.community))
    jaccard.dis.global<-rbind(jaccard.dis.global,c(jaccard.dis.temp))
  }
}

# source the function of jaccard dissimiliarity
source("JaccardDissimilarity.R")

# Calculated jaccard dissimiliarity among communities in four different environments
jaccard.dis.tropical<-Jaccard.function(130,Community.FuncType.ClimateCode) # tropical
jaccard.dis.arid<-Jaccard.function(230,Community.FuncType.ClimateCode) # arid
jaccard.dis.temperate<-Jaccard.function(330,Community.FuncType.ClimateCode) # temperate
jaccard.dis.cold<-Jaccard.function(430,Community.FuncType.ClimateCode) # cold

# pairwise t-test
t.test(jaccard.dis.tropical[,1], jaccard.dis.arid[,1])
t.test(jaccard.dis.tropical[,1], jaccard.dis.temperate[,1])
t.test(jaccard.dis.tropical[,1], jaccard.dis.cold[,1])
t.test(jaccard.dis.arid[,1], jaccard.dis.temperate[,1])
t.test(jaccard.dis.arid[,1], jaccard.dis.cold[,1])
t.test(jaccard.dis.temperate[,1], jaccard.dis.cold[,1])

# boxplot of jaccardis dissimilarity among communities in four different environments
boxplot(c(jaccard.dis.tropical,jaccard.dis.arid,jaccard.dis.temperate,jaccard.dis.cold))

# Bootstrap of Jaccard dissimilarity
source("JaccardDissimilarity_Bootstrap.R") 
boot.jaccard.dis.tropical<-Boot.Jaccard.function(130, Community.FuncType.ClimateCode, 1000)
boot.jaccard.dis.arid<-Boot.Jaccard.function(230, Community.FuncType.ClimateCode, 1000)
boot.jaccard.dis.temperate<-Boot.Jaccard.function(330, Community.FuncType.ClimateCode, 1000)
boot.jaccard.dis.cold<-Boot.Jaccard.function(430, Community.FuncType.ClimateCode, 1000)



# mean and sd of the data
mean.sd.jaccard.tropical<-cbind(mean(jaccard.dis.tropical[,1]),
                                sd(jaccard.dis.tropical[,1]),
                                length(jaccard.dis.tropical[,1]))
mean.sd.jaccard.arid<-cbind(mean(jaccard.dis.arid[,1]),
                            sd(jaccard.dis.arid[,1]),
                            length(jaccard.dis.arid[,1]))
mean.sd.jaccard.temperate<-cbind(mean(jaccard.dis.temperate[,1]),
                                 sd(jaccard.dis.temperate[,1]),
                                 length(jaccard.dis.temperate[,1]))
mean.sd.jaccard.cold<-cbind(mean(jaccard.dis.cold[,1]),
                            sd(jaccard.dis.cold[,1]),
                            length(jaccard.dis.cold[,1]))

# mean and sd of the bootstrapped data

# add all data together
mean.sd.jaccard.together<-rbind(mean.sd.jaccard.tropical,
                                mean.sd.jaccard.arid,
                                mean.sd.jaccard.temperate,
                                mean.sd.jaccard.cold)

barx <- barplot(mean.sd.jaccard.together[,1], 
                names.arg=c("Tropical","Arid","Temperate","Cold"), 
                col="grey70", 
                cex.axis=1.2,
                axis.lty=1.2, 
                cex.names=1.2,
                ylim=c(0,1),
                las=1, beside = T)
## plot the error bar
error.bar(barx, 
          mean.sd.jaccard.together[,1], 
          1.96*mean.sd.jaccard.together[,2]/sqrt(mean.sd.jaccard.together[,3]))
mtext(outer=F, side=2, line=3, text="Jaccard Dissimilarity",cex=1.5)


## calculate the differences between ecological diversity among communities (pairwise)
ed.diff.global<-data.frame()
for (q in 1:(dim(ed.global)[1]-1)) {
  for (p in q:(dim(ed.global)[1]-1)) {
    ed.diff.temp<-abs(ed.global[q,]-ed.global[p+1,])
    ed.diff.global<-rbind(ed.diff.global,c(ed.diff.temp))
  }
}


#######################################################################################################
# Barplots of ecological disparity and diversity and pair-wise jaccard dissimilarit of extant communities
#######################################################################################################
mean.ld <- t(matrix(cbind(mean(dep.tropical[,1]),
                          mean(dep.arid[,1]),
                          mean(dep.temperate[,1]),
                          mean(dep.cold[,1])),nrow=4))

sd.ld <- t(matrix(cbind(sd(dep.tropical[,1]),
                        sd(dep.arid[,1]),
                        sd(dep.temperate[,1]),
                        sd(dep.cold[,1])),nrow=4))

### Mea Diversity
mean.ld.ed <- t(matrix(cbind(mean(ed.tropical[,1]),
                             mean(ed.arid[,1]),
                             mean(ed.temperate[,1]),
                             mean(ed.cold[,1])),nrow=4))

sd.ld.ed <- t(matrix(cbind(sd(ed.tropical[,1]),
                           sd(ed.arid[,1]),
                           sd(ed.temperate[,1]),
                           sd(ed.cold[,1])),nrow=4))

# mean and sd of the jaccard dissimilarity
mean.jaccard.dis.together<-t(matrix(cbind(mean(jaccard.dis.tropical[,1]),
                                          mean(jaccard.dis.arid[,1]),
                                          mean(jaccard.dis.temperate[,1]),
                                          mean(jaccard.dis.cold[,1])), nrow=4))
sd.jaccard.dis.together<-t(matrix(cbind(sd(jaccard.dis.tropical[,1]),
                                        sd(jaccard.dis.arid[,1]),
                                        sd(jaccard.dis.temperate[,1]),
                                        sd(jaccard.dis.cold[,1])), nrow=4))

# length of each data
length.ld<-t(matrix(cbind(length(dep.tropical[,1]),
                          length(dep.arid[,1]),
                          length(dep.temperate[,1]),
                          length(dep.cold[,1])),nrow=4))

length.ld.ed<-t(matrix(cbind(length(ed.tropical[,1]),
                             length(ed.arid[,1]),
                             length(ed.temperate[,1]),
                             length(ed.cold[,1])),nrow=4))
length.jaccard.dis.together<-t(matrix(cbind(length(jaccard.dis.tropical[,1]),
                                            length(jaccard.dis.arid[,1]),
                                            length(jaccard.dis.temperate[,1]),
                                            length(jaccard.dis.cold[,1])), nrow=4))

par(mfcol=c(1,1),
    oma=c(1,1,0,1), 
    mar=c(2,4,0,1), 
    mgp=c(1,0.6,0),
    tck=-0.02)

### Without fossoil communities #### Important Codes
mean.together<-rbind(mean.ld,mean.ld.ed, mean.jaccard.dis.together*10)
sd.together<-rbind(sd.ld,sd.ld.ed, sd.jaccard.dis.together*10)
length.together<-rbind(length.ld, length.ld.ed, length.jaccard.dis.together)
barx <- barplot(mean.together, beside=T,names.arg=c("Tropical",
                                                    "Arid","Temperate","Cold"),
                ylim=c(0,17), col=c("gold3","cadetblue","darkslategrey"), cex.axis=1.2,
                axis.lty=1, cex.names=1.5,las=1, width = 0.75, xlim = c(0,13))
## plot the error bar
error.bar(barx, mean.together, 1.96*sd.together/sqrt(length.together))
mtext(outer=F, side=2, line=2, text="Mean Value",cex=1.5)
legend(0, 16,legend=c("Ecological disparity", "Ecological diversity", "Jaccard dissimilarity index"),
       fill=c("gold3","cadetblue","darkslategrey"), horiz=T,cex=1.2,box.col = "white", border="white")




##################################################
### Chinese Mesozoic mammal community analysis ###
##################################################

# Tiaojishan(TJS) community refers to Daxishan (DXS); Jiulongshan(JLS) community refers to Daohugou (DHG); based on
# Meng et al., 2015 (Arboreal docodont)
setwd("/Users/mengchen/Documents/Research Projects/MM Community Project/Database and Analysis")

ChMeMaCo<-read.csv(file="ChinaMesozoicMammalCommunity_03Feb2016.csv", header = T, sep=",")

#### Density plot of functional groups of Early Cretaceous mammal communties
par(mfcol=c(2,1),
    oma=c(3,3,1,2), 
    mar=c(2,2,1,1), 
    mgp=c(1,0.6,0),
    tck=-0.02,
    cex=0.7,
    lwd=2)

### different communities
data.dwzz<-ChMeMaCo[ChMeMaCo$Community=="DWZZ",][,4:6]   # DWZZ
write.table(data.dwzz, file="data.dwzz.csv", sep=",",row.names=F, col.names = F)
data.ljt_jsg <- ChMeMaCo[ChMeMaCo$Community=="LJT-JSG",][,4:6]   # LJT-JSG
write.table(data.ljt_jsg, file="data.ljt_jsg.csv", sep=",",row.names=F, col.names = F)
data.tjs<- ChMeMaCo[ChMeMaCo$Community=="TJS",][,4:6] # TJS
write.table(data.tjs, file="data.tjs.csv", sep=",", row.names=F, col.names = F)
data.jls<- ChMeMaCo[ChMeMaCo$Community=="JLS",][,4:6] #JLS
write.table(data.jls, file="data.jls.csv", sep=",", row.names=F, col.names = F)

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

# read data from Messel Mammalian community data
mmcd<-read.csv(file="MMCD_06042016.csv", header=F, sep=",")
eds.mmcd<-Eds.extinct(mmcd)
edv.mmcd<-Edv.extinct(mmcd)


#######################################################################################################
# Barplots of ecological disparity and diversity of both extand and extinct communities in same figure
#######################################################################################################
mean.ld <- t(matrix(cbind(mean(dep.tropical[,1]),
                          mean(dep.arid[,1]),
                          mean(dep.temperate[,1]),
                          mean(dep.cold[,1])),nrow=4))

sd.ld <- t(matrix(cbind(sd(dep.tropical[,1]),
                        sd(dep.arid[,1]),
                        sd(dep.temperate[,1]),
                        sd(dep.cold[,1])),nrow=4))

### Mea Diversity
mean.ld.ed <- t(matrix(cbind(mean(ed.tropical[,1]),
                             mean(ed.arid[,1]),
                             mean(ed.temperate[,1]),
                             mean(ed.cold[,1])),nrow=4))

sd.ld.ed <- t(matrix(cbind(sd(ed.tropical[,1]),
                           sd(ed.arid[,1]),
                           sd(ed.temperate[,1]),
                           sd(ed.cold[,1])),nrow=4))


length.ld<-t(matrix(cbind(length(dep.tropical[,1]),
                          length(dep.arid[,1]),
                          length(dep.temperate[,1]),
                          length(dep.cold[,1])),nrow=4))

length.ld.ed<-t(matrix(cbind(length(ed.tropical[,1]),
                             length(ed.arid[,1]),
                             length(ed.temperate[,1]),
                             length(ed.cold[,1])),nrow=4))
par(mfcol=c(1,2),
    oma=c(1,1,0,1), 
    mar=c(1,2,0,1), 
    mgp=c(1,0.6,0),
    tck=-0.02)

### Without fossoil communities #### Important Codes
mean.together<-rbind(mean.ld,mean.ld.ed)
sd.together<-rbind(sd.ld,sd.ld.ed)
length.together<-rbind(length.ld, length.ld.ed)
barx <- barplot(mean.together, beside=T,names.arg=c("Tropical",
                                                    "Arid","Temperate","Cold"),
                ylim=c(0,17), col=c("darkorange","steelblue3"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,10))
## plot the error bar
error.bar(barx, mean.together, 1.96*sd.together/sqrt(length.together))
mtext(outer=F, side=2, line=2, text="Mean Value",cex=1.2)
mtext(outer=F, side=3, line=-3.5, text="Extant small-bodied mammal \ncommunities", cex=1.5)
abline(h=mean(dep.global[,1]),lty=2, col="darkorange", lwd=2)
abline(h=mean(ed.global[,1]),lty=2, col="steelblue3", lwd=2)
points(1.12, 4, pch=8, col="black", cex=1.7, lwd=2,xlim=c(1,10))
points(4.13, 8.6, pch=13, col="black", cex=1.7, lwd=2.5)

###
# plots of ecological disparity and diversity of four extinct mammalian communities
mean.eds.extinct <- t(matrix(cbind(mean(eds.dwzz[,1]),
                                   mean(eds.ljt_jsg[,1]),
                                   mean(eds.tjs[,1]),
                                   mean(eds.jls[,1])),nrow=4))
sd.eds.extinct <- t(matrix(cbind(sd(eds.dwzz[,1]),
                                 sd(eds.ljt_jsg[,1]),
                                 sd(eds.tjs[,1]),
                                 sd(eds.jls[,1])),nrow=4))
# NO MEAN diversity for each extinct community; all together
all.edv.extinct<-t(matrix(cbind(edv.dwzz[,1],
                                edv.ljt_jsg[,1],
                                edv.tjs[,1],
                                edv.jls[,1])))
# length of ecological dispairity and diversity
length.eds.extinct<-t(matrix(cbind(length(eds.dwzz[,1]),
                                   length(eds.ljt_jsg[,1]),
                                   length(eds.tjs[,1]),
                                   length(eds.jls[,1])),nrow=4))

mean.eds.edv.extinct<-rbind(mean.eds.extinct, all.edv.extinct)
sd.eds.edv.extinct<-rbind(sd.eds.extinct,c(0,0,0,0))
length.eds.edv.extinct<-rbind(length.eds.extinct, c(1,1,1,1))
barx <- barplot(mean.eds.edv.extinct, beside=T,names.arg=c("DWZZ","LJU-JSG","TJS","JLS"),
                ylim=c(0,17), col=c("darkorange","steelblue3"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,10))
## plot the error bar
error.bar(barx, mean.eds.edv.extinct, 1.96*sd.eds.edv.extinct/sqrt(length.eds.edv.extinct))
abline(h=mean(dep.global[,1]),lty=2, col="darkorange", lwd=2)
abline(h=mean(ed.global[,1]),lty=2, col="steelblue3", lwd=2)
legend(4, 13.5,legend=c("Ecological disparity", "Ecological diversity"),
       fill=c("darkorange","steelblue3"),horiz=F,cex=1.2,
       box.col = "white",border="white")
mtext(outer=F, side=3, line=-3, text="Mesozoic mammal communities", cex=1.5)
text(x=10, y=mean(dep.global[,1])+0.5, mean(dep.global[,1]),col="darkorange")
text(x=10, y=mean(ed.global[,1])+0.5, mean(ed.global[,1]),col="steelblue3")

##################################################
# Jaccard dissimiliarity among extinct communities
##################################################
Extinct.Community.FuncType<-read.csv(file="ExtinctCommunityEdData_03Feb2016.csv", header = T, sep=",")
Extinct.community<-unique(Extinct.Community.FuncType$CommunityNo)
jaccard.dis<-data.frame()

# jaccard function
Extinct.Community.Jaccard.dis<- function (input.data) {
  for (g in 1:(length(Extinct.community)-1)) {
    for (h in g:(length(Extinct.community)-1)) {
      func.community.temp1<-unique(input.data[input.data$CommunityNo==g,]$FuncType)
      func.community.temp2<-unique(input.data[input.data$CommunityNo==h+1,]$FuncType)
      # jaccarad disimiliarity equation
      jaccard.dis.temp<-1-length(intersect(func.community.temp1,func.community.temp2))/length(union(func.community.temp1,func.community.temp2))
      jaccard.dis<-rbind(jaccard.dis,c(jaccard.dis.temp))
    }
  }
  return(jaccard.dis)
}

# Jaccard dissimilarity among four extinct communities
Extinct.Community.Jaccard.dis(Extinct.Community.FuncType)

# Jaccard dissimilarity between extinct communities and extant communities from four environments

# jaccard function
Community.Jaccard.dis <- function (extinct.data, extant.data, code) {
  jaccard.dis<-data.frame()
  if (code < 200) {
    extant.data<-extant.data[extant.data$ClimateCode < 200,] # tropical
  } 
  if (code > 200 && code < 300) {
    extant.data<-extant.data[(extant.data$ClimateCode >200) & (extant.data$ClimateCode < 300),] # arid
  } 
  if (code > 300 && code < 400) {
    extant.data<-extant.data[(extant.data$ClimateCode >300) & (extant.data$ClimateCode < 400),]  # temperate
  } 
  if (code > 400) {
    extant.data<-extant.data[extant.data$ClimateCode > 400,] # cold
  }
  
  for (g in unique(extinct.data[,1])) {
    for (h in unique(extant.data[,1])) {
      func.community.temp1<-unique(extinct.data[extinct.data$CommunityNo==g,]$FuncType) # extinct community
      func.community.temp2<-unique(extant.data[extant.data$No==h,]$FuncType) # extant community
      # jaccarad disimiliarity equation
      jaccard.dis.temp<-1-length(intersect(func.community.temp1,func.community.temp2))/
        length(union(func.community.temp1,func.community.temp2))
      jaccard.dis<-rbind(jaccard.dis,c(jaccard.dis.temp))
    }
  }
  return(jaccard.dis)
}

# Extinct commnities vs tropical extant commnities
tropical.community<-Community.FuncType.ClimateCode[Community.FuncType.ClimateCode$ClimateCode<200,]$No
tropical.community<-unique(tropical.community)
Extinct.tropical.jaccard.dis<-Community.Jaccard.dis(Extinct.Community.FuncType, Community.FuncType.ClimateCode ,130)
Extinct.tropical.jaccard.dis<-matrix(as.matrix(Extinct.tropical.jaccard.dis), nrow=15, ncol=5, byrow = F)
write.table(Extinct.tropical.jaccard.dis, file = "Extinct.tropical.jaccard.dis.csv", sep=",",
            row.names=tropical.community, col.names=c("DWZZ","LJT-JSG","TJS","JLS", "MSL"))
# Extinct communitis vs aird extant communities
arid.community<-Community.FuncType.ClimateCode[Community.FuncType.ClimateCode$ClimateCode>200 & 
                                                     Community.FuncType.ClimateCode$ClimateCode<300,]$No
arid.community<-unique(arid.community)
Extinct.arid.jaccard.dis<-Community.Jaccard.dis(Extinct.Community.FuncType, Community.FuncType.ClimateCode ,230)
Extinct.arid.jaccard.dis<-matrix(as.matrix(Extinct.arid.jaccard.dis), nrow=10, ncol=5, byrow = F) 
write.table(Extinct.arid.jaccard.dis, file = "Extinct.arid.jaccard.dis.csv", sep=",", 
            row.names=arid.community, col.names=c("DWZZ","LJT-JSG","TJS","JLS", "MSL"))
# Extinct communitis vs temperate extant communities
temperate.community<-Community.FuncType.ClimateCode[Community.FuncType.ClimateCode$ClimateCode>300 & 
                                                 Community.FuncType.ClimateCode$ClimateCode<400,]$No
temperate.community<-unique(temperate.community)
Extinct.temperate.jaccard.dis<-Community.Jaccard.dis(Extinct.Community.FuncType, Community.FuncType.ClimateCode ,330)
Extinct.temperate.jaccard.dis<-matrix(as.matrix(Extinct.temperate.jaccard.dis), nrow=15, ncol=5, byrow = F)
write.table(Extinct.temperate.jaccard.dis, file = "Extinct.temperate.jaccard.dis.csv", sep=",", 
            row.names=temperate.community, col.names=c("DWZZ","LJT-JSG","TJS","JLS", "MSL"))
# Extinct communitis vs cold extant communities
cold.community<-Community.FuncType.ClimateCode[Community.FuncType.ClimateCode$ClimateCode>400,]$No
cold.community<-unique(cold.community)
Extinct.cold.jaccard.dis<-Community.Jaccard.dis(Extinct.Community.FuncType, Community.FuncType.ClimateCode ,430)
Extinct.cold.jaccard.dis<-matrix(as.matrix(Extinct.cold.jaccard.dis), nrow=5, ncol=5, byrow = F)
write.table(Extinct.cold.jaccard.dis, file = "Extinct.cold.jaccard.dis.csv", sep=",", 
            row.names=cold.community, col.names=c("DWZZ","LJT-JSG","TJS","JLS","MSL"))


# DWZZ
dwzz.mean.jaccard.dis<-cbind(mean(Extinct.tropical.jaccard.dis[,1]), mean(Extinct.arid.jaccard.dis[,1]),
                             mean(Extinct.temperate.jaccard.dis[,1]), mean(Extinct.cold.jaccard.dis[,1]))
dwzz.sd.jaccard.dis<-cbind(sd(Extinct.tropical.jaccard.dis[,1]), sd(Extinct.arid.jaccard.dis[,1]),
                           sd(Extinct.temperate.jaccard.dis[,1]), sd(Extinct.cold.jaccard.dis[,1]))
# ljt_jsg
ljt_jsg.mean.jaccard.dis<-cbind(mean(Extinct.tropical.jaccard.dis[,2]), mean(Extinct.arid.jaccard.dis[,2]),
                                mean(Extinct.temperate.jaccard.dis[,2]), mean(Extinct.cold.jaccard.dis[,2]))
ljt_jsg.sd.jaccard.dis<-cbind(sd(Extinct.tropical.jaccard.dis[,2]), sd(Extinct.arid.jaccard.dis[,2]),
                              sd(Extinct.temperate.jaccard.dis[,2]), sd(Extinct.cold.jaccard.dis[,2]))
# tjs
tjs.mean.jaccard.dis<-cbind(mean(Extinct.tropical.jaccard.dis[,3]), mean(Extinct.arid.jaccard.dis[,3]),
                            mean(Extinct.temperate.jaccard.dis[,3]), mean(Extinct.cold.jaccard.dis[,3]))
tjs.sd.jaccard.dis<-cbind(sd(Extinct.tropical.jaccard.dis[,3]), sd(Extinct.arid.jaccard.dis[,3]),
                          sd(Extinct.temperate.jaccard.dis[,3]), sd(Extinct.cold.jaccard.dis[,3]))
# jls
jls.mean.jaccard.dis<-cbind(mean(Extinct.tropical.jaccard.dis[,4]), mean(Extinct.arid.jaccard.dis[,4]),
                            mean(Extinct.temperate.jaccard.dis[,4]), mean(Extinct.cold.jaccard.dis[,4]))
jls.sd.jaccard.dis<-cbind(sd(Extinct.tropical.jaccard.dis[,4]), sd(Extinct.arid.jaccard.dis[,4]),
                          sd(Extinct.temperate.jaccard.dis[,4]), sd(Extinct.cold.jaccard.dis[,4]))

#### Plots the jaccard dissimiliary index
layout(matrix(c(1,2,1,3,1,4,1,5),4,2, 
              byrow = T), 
       width=c(1,1.5))
# DWZZ
low.diff.jaccard.dis<-dwzz.mean.jaccard.dis-dwzz.sd.jaccard.dis
upper.diff.jaccard.dis<-dwzz.mean.jaccard.dis+dwzz.sd.jaccard.dis
plot(dwzz.mean.jaccard.dis[1,],col="darkorange3", type="b", pch=19, lty=2, lwd=1.5, 
     ylim=c(0.02,1), xaxt="none", yaxt="none",xlab="",ylab="",cex=2)
arrows(1:4,low.diff.jaccard.dis, 1:4, upper.diff.jaccard.dis, length=0.05,angle = 90, code=3, col="darkorange3", lwd=1.5)
# 
low.diff.jaccard.dis<-ljt_jsg.mean.jaccard.dis-ljt_jsg.sd.jaccard.dis-0.25
upper.diff.jaccard.dis<-ljt_jsg.mean.jaccard.dis+ljt_jsg.sd.jaccard.dis-0.25
points(ljt_jsg.mean.jaccard.dis[1,]-0.25,col="grey30", type="b", pch=19, lty=2, lwd=1.5, cex=2)
arrows(1:4,low.diff.jaccard.dis, 1:4, upper.diff.jaccard.dis, length=0.05,angle = 90, code=3, col="grey30", lwd=1.5)

low.diff.jaccard.dis<-tjs.mean.jaccard.dis-tjs.sd.jaccard.dis-0.5
upper.diff.jaccard.dis<-tjs.mean.jaccard.dis+tjs.sd.jaccard.dis-0.5
points(tjs.mean.jaccard.dis[1,]-0.5,col="green3", type="b", pch=19, lty=2, lwd=1.5, cex=2)
arrows(1:4,low.diff.jaccard.dis, 1:4, upper.diff.jaccard.dis, length=0.05,angle = 90, code=3, col="green3", lwd=1.5)

low.diff.jaccard.dis<-jls.mean.jaccard.dis-jls.sd.jaccard.dis-0.75
upper.diff.jaccard.dis<-jls.mean.jaccard.dis+jls.sd.jaccard.dis-0.75
points(jls.mean.jaccard.dis[1,]-0.75,col="dodgerblue3", type="b", pch=19, lty=2, lwd=1.5, cex=2)
arrows(1:4,low.diff.jaccard.dis, 1:4, upper.diff.jaccard.dis, length=0.05,angle = 90, code=3, col="dodgerblue3", lwd=1.5)

abline(h=0.25,lty=2, col="grey40")
abline(h=0.5,lty=2, col="grey40")
abline(h=0.75,lty=2, col="grey40")
abline(h=1, lty=2, col="grey40")
abline(h=0.025)
abline(h=.275)
abline(h=0.525)
abline(h=0.775)
mtext(outer = F, side=1, line=2, text="Tropical        Arid       Temperate        Cold", cex=1.2)

# DWZZ
barx <- barplot(mean.eds.edv.extinct[,1], beside=T,
                ylim=c(0,9), col=c("gold3","cadetblue"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,15))
## plot the error bar
error.bar(barx, mean.eds.edv.extinct[,1], 1.96*sd.eds.edv.extinct[,1]/sqrt(length.eds.edv.extinct[,1]))
# LJT-JSG
barx <- barplot(mean.eds.edv.extinct[,2], beside=T,
                ylim=c(0,9), col=c("gold3","cadetblue"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,15))
## plot the error bar
error.bar(barx, mean.eds.edv.extinct[,2], 1.96*sd.eds.edv.extinct[,2]/sqrt(length.eds.edv.extinct[,2]))
#TJS
barx <- barplot(mean.eds.edv.extinct[,3], beside=T,
                ylim=c(0,9), col=c("gold3","cadetblue"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,15))
## plot the error bar
error.bar(barx, mean.eds.edv.extinct[,3], 1.96*sd.eds.edv.extinct[,3]/sqrt(length.eds.edv.extinct[,3]))
#JLS
barx <- barplot(mean.eds.edv.extinct[,4], beside=T,
                ylim=c(0,9), col=c("gold3","cadetblue"), cex.axis=1.2,
                axis.lty=1, cex.names=1.2,las=1, width = 0.75, xlim = c(0,15))
## plot the error bar
error.bar(barx, mean.eds.edv.extinct[,4], 1.96*sd.eds.edv.extinct[,4]/sqrt(length.eds.edv.extinct[,4]))


#########################################################
# When two communities have the joint absencees of same functional types,
# this might mean that the communities are ecologicall similar. On the other hand. if a functional
# type has a highly clumped distribution, or is simply rare, then join basences might arise through chance and say
# nothing abou the suitability of a given community for a species, the similarity among the community needs of function
# type or the cological similarity of communities
# this is similart to the statement in Zuur et al., 2009: page 8, left bottom

# Function type data table
FuncType.occur<-xtabs(~Community.FuncType$No+Community.FuncType$FuncType)

# get rid of any values that might be greater than one

for(i in 1:nrow(FuncType.occur)){
  for(j in 1:ncol(FuncType.occur)){
    if(FuncType.occur[i,j]>1){
      FuncType.occur[i,j]<-1
    }
  }
}

library(corrgram) # this is based on the PC1 and PC2 scores; doest fit for the data I have here.
corrgram(FuncType.occur[,1:107], order = T)

# rest of the below doesnt work
FuncType.occur.data<-FuncType.occur[,1:107]

A <- matrix(nrow = N, ncol = N)

for (i in 1:N){
  for (j in 1:N){
    A[i,j] <- sum(FuncType.occur.data[i]==0  & FuncType.occur.data[j]==0, na.rm=TRUE)
  }}


A1 <- A/2035
print(A1, digits = 2)
rownames(A1) <- AllNames
colnames(A1) <- AllNames


library(lattice)

panel.corrgram.2 <- function(x, y, z, subscripts, at = pretty(z), scale = 0.8, ...)
{
  require("grid", quietly = TRUE)
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  zcol <- level.colors(z, at = at, ...)
  for (i in seq(along = z))
  {
    lims <- range(0, z[i])
    tval <- 2 * base::pi *
      seq(from = lims[1], to = lims[2], by = 0.01)
    grid.polygon(x = x[i] + .5 * scale * c(0, sin(tval)),
                 y = y[i] + .5 * scale * c(0, cos(tval)),
                 default.units = "native",
                 gp = gpar(fill = zcol[i]))
    grid.circle(x = x[i], y = y[i], r = .5 * scale,
                default.units = "native")
  }
}




levelplot(A1,xlab=NULL,ylab=NULL,
          at=do.breaks(c(0.5,1.01),101),
          panel=panel.corrgram.2,
          scales=list(x=list(rot=90)),
          colorkey=list(space="top"),
          col.regions=colorRampPalette(c("red","white","blue")))


#Grey colours
levelplot(A1,xlab=NULL,ylab=NULL,
          at=do.breaks(c(0.5,1.01),101),
          panel=panel.corrgram.2,
          scales=list(x=list(rot=90)),
          colorkey=list(space="top"),
          col.regions=colorRampPalette(c(grey(0.8),grey(0.5),grey(0.2))))

