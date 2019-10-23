# Geographic plot codes written by Meng Chen during the doctoral dissertation.

# 22April,2014
# revised on 19Jan2016
# reorganized on 23Oct2019

###########################################################
# Geographic plots of the small-bodied mammal communities #
# around the world.                                       #
# Author: Meng Chen                     Date: 12Jan2016   #
###########################################################

# require packages for mapping available
library(maps)
library(mapdata)
library(mapproj)
library(vegan)
library(MASS)
library(ggplot2)
library(ggmap)
library(ggthems)

setwd("your directory")

### the number of species in each commmunity
mcm.data <- read.csv(file="small mammal community data", sep=",", header=T)
mcm.data <- mcm.data[-c(546:573, 821:823, 830:831),]
n.species.each.community <- table(mcm.data$Community.No)
n.species.each.community <- as.data.frame(n.species.each.community)
colnames(n.species.each.community) <- c("Community", "No_species")

# GPS coordinates of small mammal communities
coords <- read.csv(file="GPS data", header=TRUE, sep = ",")
coords <- cbind(coords, n.species.each.community )
coords <- as.data.frame(coords)

# geographic plots using gglot2
# basic maps
world <- map_data("world")
wrold <- world[world$region != "Antarctica",]
gg <- ggplot()

# 2 habitats
gg.habitat <- gg + geom_map(data=world, map = world, aes(x=long, y=lat, map_id=region), color="black", fill="grey60", size=0.1, alpha=0.5) +
  geom_point(data = coords, aes(x=Lon, y=Lat, color=Habitat), size=3, alpha=1) +
  scale_color_manual(values=c("grey60", "black")) +
  theme(legend.position="bottom") +
  theme_bw()
ggsave("habitats.pdf", plot=last_plot(), width=9, height=5, scale=1, dpi=600)

# 4 climates ---- used in the manuscript Fig. 1
gg.climate <- gg + geom_map(data=world, map = world, aes(x=long, y=lat, map_id=region), color="black", fill="grey60", size=0.1, alpha=0.5) +
  ggplot(coords, aes(x=Lon, y=Lat)) +
  geom_point(aes(fill=Climate, color=Climate, size=No_species), alpha=0.8) +
  scale_fill_manual(values=c("#584F59","dodgerblue3","green3","darkorange","grey30")) +
  scale_color_manual (values=c("grey80","steelblue1","palegreen","orange","grey30")) +
  theme(legend.position="bottom") +
  theme_bw()
ggsave("climates.pdf", plot=last_plot(), width=9, height=5, scale=1, dpi=600)

# 8 vegetations
gg.vegetation <- gg + geom_map(data=world, map = world, aes(x=long, y=lat, map_id=region), color="black", fill="grey60", size=0.1, alpha=0.5) +
   geom_point(data = coords, aes(x=Lon, y=Lat, color=Vegetation), size=3, alpha=1) +
   scale_color_manual(values=c("darkviolet","khaki","goldenrod", "firebrick1", "darkgoldenrod1", "green3", "red", "orangered")) +
   theme(legend.position="bottom") +
   theme_bw()
ggsave("vegetations.pdf", plot=last_plot(), width=9, height=5, scale=1, dpi=600)


