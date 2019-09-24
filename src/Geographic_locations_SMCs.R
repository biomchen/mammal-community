## The code is written by Meng Chen
## for dissertation project III: paleoecology
## of Early Cretaceous mammals
## 22April,2014

## revised on 19Jan2016

## require packages for mapping available 
## mammal communities around the world
library(maps)       #basic mapping functions and some data
library(mapdata)    #some additional hires data
library(mapproj)
library(vegan)
library(MASS)
library(ggplot2)
library(ggmap)
library(ggthems)

setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis")

### the number of species in each commmunity
mcm.habitat<-read.csv(file="EcoHabitatAll_6June2017.csv", sep=",", header=T)
mcm.habitat<- mcm.habitat[-c(546:573, 821:823, 830:831),]
n.species.each.community <- table(mcm.habitat$Community.No)
n.species.each.community <- as.data.frame(n.species.each.community)
colnames(n.species.each.community) <- c("Community", "No_species")


setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis/Geographic Plots")

### GPS coordinates of each community
coords<-read.csv(file="GPS_coordinate_13July2017.csv", header=TRUE, sep = ",")
coords<-cbind(coords,n.species.each.community )
coords<-as.data.frame(coords)

### gglot

world <- map_data("world")
wrold <- world[world$region != "Antarctica",]
gg<- ggplot()
## habitat
gg.habitat <- gg + geom_map(data=world, map = world, aes(x=long, y=lat, map_id=region), color="black", fill="grey60", size =0.1, alpha = 0.5)
gg.habitat <- gg.habitat + geom_point(data = coords, aes(x=Lon, y=Lat, color=Habitat), size=3, alpha=1)
gg.habitat <- gg.habitat + scale_color_manual(values = c("grey60", "black"))
## climate

gg.climate <- gg + geom_map(data=world, map = world, aes(x=long, y=lat, map_id=region), color="black", fill="grey60", size =0.1, alpha = 0.5) +
  ggplot(coords, aes(x=Lon, y=Lat)) +
  geom_point(aes(fill=Climate, color=Climate, size=No_species), alpha=0.8) +
  scale_fill_manual(values = c("#584F59","dodgerblue3","green3","darkorange","grey30")) +
  scale_color_manual (values = c("grey80","steelblue1","palegreen","orange","grey30")) +
  theme(legend.position = "bottom") + 
  theme_bw()

ggsave("gg.climate.pdf", plot = last_plot(), width = 9, height = 5, scale = 1, dpi = 600)

## vegetation
gg.vegetation <- gg + geom_map(data=world, map = world, aes(x=long, y=lat, map_id=region), color="black", fill="grey60", size =0.1, alpha = 0.5)
gg.vegetation <- gg.vegetation + geom_point(data = coords, aes(x=Lon, y=Lat, color=Vegetation), size=3, alpha=1)
gg.vegetation <- gg.vegetation + scale_color_manual(values = c("darkviolet","khaki","goldenrod", "firebrick1", "darkgoldenrod1", "green3","red", "orangered"))


##############################
## ## google map plots
## I dont know why is offset between the map and points
world.map <- qmap("world",  zoom =1)
world.map + geom_point(data = coords, mapping = aes(x=Lon, y=Lat, color=Habitat), size=3, alpha=1)
world.map + geom_point(data = coords, mapping = aes(x=Lon, y=Lat, color=Climate), size=3, alpha=1)
world.map + geom_point(data = coords, mapping = aes(x=Lon, y=Lat, color=Vegetation), size=3, alpha=1)

####################################
## modifed original plots
## habitat color
col.list.habitat <- c("grey60", "black")
palette(col.list.habitat)
pairs(coords[,5:6], col=coords[,2])
coords$Habitat.color <- palette()[coords[,2]]
## check if it is right
table(coords[,2], palette()[coords[,2]])

## climate color
col.list.climate <- c("grey20","dodgerblue3","green3","darkorange")
palette(col.list.climate)
pairs(coords[,5:6], col=coords[,3])
coords$Climate.color <- palette()[coords[,3]]
## check if it is right
table(coords[,3], palette()[coords[,3]])

## vegetation
col.list.vegetation <- c("darkviolet", "khaki","goldenrod","firebrick1","darkgoldenrod1","green3", "red", "orangered")
palette(col.list.vegetation)
pairs(coords[,5:6], col=coords[,4])
coords$Vegetation.color <- palette()[coords[,4]]
## check if it is right
table(coords[,4], palette()[coords[,4]])

#### original plot
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(0,0,0,0))

### habitat
map('worldHires',col="grey30", fill = TRUE)
lat<-c(-60, -40, -20, 0, 20, 40, 60, 80)

for (j in 1:8) {
  abline(h=lat[j], lty=5, col="grey60")
}
points(coords$Lon,coords$Lat,
			 pch=21,bg=coords$Habitat.color, cex=coords$No_species/max(coords$No_species)*5)

### climate
map('worldHires',col="grey30", fill = TRUE)
lat<-c(-60, -40, -20, 0, 20, 40, 60, 80)

for (j in 1:8) {
  abline(h=lat[j], lty=5, col="grey60")
}
points(coords$Lon,coords$Lat,
       pch=21,bg=coords$Climate.color, cex=coords$No_species/max(coords$No_species)*5)

### vegetation
map('worldHires',col="grey30", fill = TRUE)
lat<-c(-60, -40, -20, 0, 20, 40, 60, 80)

for (j in 1:8) {
  abline(h=lat[j], lty=5, col="grey60")
}
points(coords$Lon,coords$Lat,
       pch=21,bg=coords$Vegetation.color, cex=coords$No_species/max(coords$No_species)*5)


points(x=c(-145,-145,-145,-145,
           -145,-145,-145,-145),
       y=c(8,4,-1,-8,-16,-26,-38,-51),pch=21,bg="white",
       cex=c(3/30*4,7/30*5,11/30*5,
             15/30*5,19/30*5,23/30*5,
             27/30*5,30/30*5))

text(x=c(-145,-145,-145,-145,-145,
         -145,-145,-145,-145),
     y=c(14,8,4,-1,-8,-16,-26,-38,-51),
     labels=c("No. of Species","3","7","11",
              "15","19","23","27","30"),
     offset=1,cex=0.6,pos=4)

mtext(outer=T, side=1, line=1, 
      text="45 Worldwide-Sampled Small Mammal Commuity",
      cex=1,font=1,las=1)


###################################################
## Plot the map of world climte from 1977 to 2000
require(rworldmap)
require(maptools)
require(sp)
require(fields)
World_Climate_1977_2000 <- 'Koeppen-geiger-ASCII.txt'
#read in data which is as lon,lat,catID 
dataframe<-read.table(World_Climate_1977_2000,header=TRUE,as.is=TRUE) 
#convert to sp SpatialPointsDataFrame 
coordinates(dataframe) = c("Lon", "Lat")
# promote to SpatialPixelsDataFrame 
gridded(dataframe) <- TRUE
# promote to SpatialGridDataFrame
sGDF = as(dataframe, "SpatialGridDataFrame")
col<-read.csv(file="Color for climate.csv",as.is=TRUE)

#plotting map
mapDevice() 
#create world map shaped window 
mapParams <- mapGriddedData(sGDF,
														catMethod='categorical',
														addBorders = "none",
														addLegend=FALSE,
														colourPalette=col$ColorCode)
#adding formatted legend
do.call(addMapLegendBoxes,
				c(mapParams,
					cex=0.8,
					ncol=11,
					x='bottom',
					title='Koeppen-Geiger Climate Zones'))

quartz.save("Koeppen-Geiger Classification.pdf","pdf") 


max<-max(coords$NoSpecies)

points(coords$Lon,coords$Lat,
			 pch=21, bg=coords$Color,cex=coords$NoSpecies/max*4)
# text(coords$Longtitude,coords$Latitude,pos = 2, cex=0.5)

points(x=c(-145,-145,-145,-145,
					 -145,-145,-145,-145),
			 y=c(8,4,-1,-8,-16,-26,-38,-51),pch=21,bg="white",
			 cex=c(3/30*4,7/30*4,11/30*4,
			 			15/30*4,19/30*4,23/30*4,
			 			27/30*4,30/30*4))

text(x=c(-145,-145,-145,-145,-145,
         -145,-145,-145,-145),
     y=c(14,8,4,-1,-8,-16,-26,-38,-51),
     labels=c("No. of Species","3","7","11",
              "15","19","23","27","30"),
     offset=1,cex=0.6,pos=4)

