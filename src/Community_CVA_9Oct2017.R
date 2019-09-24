###04/17/2013

##DFA on osteological indices(<0.001)
library(vegan)
library(MASS)
library(scatterplot3d)
library(permute)
library(calibrate)
library(car)
library(grid)
library(vcd)
library(rgl)
library(ggplot2)
library(gridExtra)
library(cowplot)

source("biostats.R")  ##helps with data screening, compiled by K. McGarigal at UMass

setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis")

## Reading the MasterSheet data
Community.proportion.data <- read.csv("EcoHabitat_Proportion_MasterSheet_6June2017.csv", header =T, sep=",")
## No of extant small-bodied mammal communities
n.community<-dim(Community.proportion.data)[1]

## Habitat data
Community.proportion.habitat <- Community.proportion.data[order(Community.proportion.data$Habitat),]
## Climate data
Community.proportion.climate <- Community.proportion.data[order(Community.proportion.data$Climate),]
## Vegetation data
Community.proportion.vegetation <- Community.proportion.data[order(Community.proportion.data$Vegetation),]

# extinct mammalian community data
ChMeMaCo<-read.csv(file="ChinaMesozoicMammalCommunity_28June2018.csv", header = T, sep=",")

### proportion data of extinct mamamlian communities

extinct.data.proportion <- cbind(table(ChMeMaCo[,c(1,4)])/rowSums(table(ChMeMaCo[,c(1,4)])), 
                                 table(ChMeMaCo[,c(1,5)])/rowSums(table(ChMeMaCo[,c(1,5)])), 
                                 table(ChMeMaCo[,c(1,6)])/rowSums(table(ChMeMaCo[,c(1,6)])))

colnames(extinct.data.proportion) <- c("BS1", "BS2", "BS3", "BS4", "BS5",
                                       "DP1", "DP2", "DP3", "DP6",
                                       "LM1", "LM2", "LM3", "LM4", "LM5", "LM6", "LM7", "LM8")

### non-edited proportion data of extinct mammalian communities
#write.csv(extinct.data.proportion, file = "Extinct_Mammal_Proportion_NonPrep_9July2018.csv")
### reading the edited extinct mammalian communities data
extinct.data<-read.csv("Extinct_Mammal_Proportion_9July2018.csv", header = T, sep=",")

########################
## OPEN vs Close habitat
########################
## combine data
combined.data.habitat<-rbind(Community.proportion.habitat[,c(2,5:23)], extinct.data[,c(2,5:23)])

# number of close habitat
n.close<-length(which(Community.proportion.habitat$Habitat=="Close"))
# number of open habitat
n.open<-length(which(Community.proportion.habitat$Habitat=="Open"))

## Discriminant Functional Analysis on Significant Indices(<0.001)
sigdata.habitat<-combined.data.habitat

# Climate<-combined.data[,1]
Habitat<-combined.data.habitat[,1]

# ## data transformation transformed data = arcsin (sqrt(original data))
sigdata.habitat<-cbind(Habitat, asin(sqrt(sigdata.habitat[,2:20])))

## apply X 100 to the precentage data
# sigdata<-cbind(Habitat, sigdata[,2:20]*100)

Strudata.habitat<-sigdata.habitat[1:n.community,2:20]

## Testing for significant differences between groups
sigindex.habitat <- as.matrix(sigdata.habitat[,c(2:20)])
sigoutput.habitat <- manova(sigindex.habitat~sigdata.habitat$Habitat)
sigsum.wilks.habitat<-summary(sigoutput.habitat, test="Wilks")
sigsum.avo.habitat<-summary.aov(sigoutput.habitat)

## DFA based on habitat types
siglda.mod.habitat <- lda(Habitat ~ ., data = sigdata.habitat[1:n.community,], prior=c(1/2, 1/2))
siglda.res.habitat <- predict(siglda.mod.habitat, sigdata.habitat)

### plots the DF1 against the latitudes
DF1_two_habitats <- cbind(Community.proportion.habitat[,1:2], 
                          siglda.res.habitat$x[-(99:103),])
DF1_two_habitats <- DF1_two_habitats[order(DF1_two_habitats$Community),]
colnames(DF1_two_habitats) <- c("Community", "Habitat", "DF1")
write.csv(DF1_two_habitats, file = "DF1_two_habitat_30June2017.csv")

## eigenvaules
siglda.mod.habitat$svd

## First two coefficient of the LDA of the osteological indices

table_2habitat<-table(sigdata.habitat$Habitat[1:n.community], siglda.res.habitat$class[1:n.community])
write.table(table_2habitat, file = "table_2habitat.csv", sep = ",")

cols.habitat=c(rep("grey30", n.close),  # close
               rep("black", n.open),  # open
               rep("grey30", 5))           


sybs.habitat=c(rep(4, n.close),       # close   
               rep(15, n.open),         # open     
               rep(16,5))    


##############################
##  Four climatic types
##############################
combined.data.climate<-rbind(Community.proportion.climate[,c(3,5:23)], extinct.data[,c(3,5:23)])

# number of each climate
n.tropical<-length(which(Community.proportion.climate$Climate=="Tropical"))
n.arid<-length(which(Community.proportion.climate$Climate=="Arid"))
n.temperate<-length(which(Community.proportion.climate$Climate=="Temperate"))
n.cold<-length(which(Community.proportion.climate$Climate=="Cold"))

## Discriminant Functional Analysis on Significant Indices(<0.001)
sigdata.climate<-combined.data.climate

Climate<-combined.data.climate[,1]

# ## data transformation transformed data = arcsin (sqrt(original data))
sigdata.climate<-cbind(Climate, asin(sqrt(sigdata.climate[,2:20])))

## apply X 100 to the precentage data
# sigdata<-cbind(Climate, sigdata[,2:20]*100)

Strudata.climate<-sigdata.climate[1:n.community,2:20]

## Testing for significant differences between groups
sigindex.climate <- as.matrix(sigdata.climate[,c(2:20)])
sigoutput.climate <- manova(sigindex.climate~sigdata.climate$Climate)
sigsum.wilks.climate<-summary(sigoutput.climate, test="Wilks")
sigsum.avo.climate<-summary.aov(sigoutput.climate)

## Discriminant functional analyses
siglda.mod.climate <- lda(Climate ~ ., data = sigdata.climate[1:n.community,], prior=c(1/4, 1/4, 1/4, 1/4))
siglda.res.climate <- predict(siglda.mod.climate, sigdata.climate)

## eigenvaules
siglda.mod.climate$svd

## First two coefficient of the LDA of the osteological indices
coemtr.climate<-cbind(siglda.mod.climate$scaling[,1],siglda.mod.climate$scaling[,2])

table_4environments<-table(sigdata.climate$Climate[1:n.community], siglda.res.climate$class[1:n.community])
write.table(table_4environments, file = "table_4environments.csv", sep = ",")

cols.climate=c(rep("white",n.arid),      ## Arid
               rep("white",n.cold),         ## Cold
               rep("white",n.temperate),         ## Temperate
               rep("white",n.tropical),      ## Tropical
               rep("red", 5))

bg.climate <-c(rep("#59505A",n.arid),      ## Arid
               rep("dodgerblue3",n.cold),         ## Cold
               rep("green3",n.temperate),         ## Temperate
               rep("darkorange",n.tropical),      ## Tropical
               rep("red", 5))           

sybs.climate=c(rep(23, n.arid),
               rep(22, n.cold),              
               rep(24, n.temperate),          
               rep(21, n.tropical),
               rep(4,5))


## Structual Matrix of DFA 
StruCol.climate<-lda.structure(siglda.res.climate$x[1:n.community,],Strudata.climate) 
ld1.climate<-as.matrix(StruCol.climate$LD1)
ld2.climate<-as.matrix(StruCol.climate$LD2)
ld3.climate<-as.matrix(StruCol.climate$LD3)

write.csv(StruCol.climate, file = "StructMatrix_four_environ.csv")
S.climate<-as.matrix(cbind(ld1.climate,ld2.climate,ld3.climate))

###
par(mfrow=c(1,1),oma=c(3,3,0,0),mar=c(2,2,0,1),
    mgp=c(3,0.5,0),tck=-0.005,cex=0.8,cex.axis=0.7)

## Plot DF1 vs DF2
plot(siglda.res.climate$x[,1],
     siglda.res.climate$x[,2],
     ann=F,
     pch=sybs.climate,
     col=cols.climate,
     bg=bg.climate,
     las=1,
     frame=F,cex=2.5)

abline(v=0,lty=3)
abline(h=0,lty=3)
mtext(outer=F, side=2, line=2, text="DF2 (31.75%)",cex=1)
textxy(siglda.res.climate$x[,1], siglda.res.climate$x[,2], cex=1, labs = Community.proportion.climate$Community)

##Structure matrix plot CF1 and CF2
plot(S.climate[,1:2], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:19) {
  arrows(0, 0, as.numeric(S.climate[i,1]), as.numeric(S.climate[i,2]), lty=1, lwd=1, 
         col="chocolate1", length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S.climate[,1]), as.numeric(S.climate[,2]), cex=1, labs = rownames(S.climate))
mtext(outer=F, side=3, line=2, text="Coefficient",las=1,cex=1)

## Plot CF1 vs CF3
plot(siglda.res.climate$x[,1],
     siglda.res.climate$x[,3],
     ann=F,
     pch=sybs.climate,
     col=cols.climate,
     bg=bg.climate,
     las=1,
     frame=F,cex=2.5)

mtext(outer=F, side=1, line=1, text="DF1 (58.19%)",cex=1)
mtext(outer=F, side=2, line=2, text="DF3 (10.89%)",cex=1)
textxy(siglda.res.climate$x[,1], siglda.res.climate$x[,3], cex=1, labs = Community.proportion.climate$Community)

##Structure matrix plot CF1 and CF3
plot(S.climate[,c(1,3)], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:19){
  arrows(0,0,as.numeric(S.climate[i,1]),as.numeric(S.climate[i,3]), 
         lty=1, 
         lwd=1, 
         col="chocolate1",
         length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S.climate[,1]), as.numeric(S.climate[,3]),cex=1, labs=rownames(S.climate))
mtext(outer=F, side=1, line=2, text="DF1",las=1,cex=1)

## Plot DF2 vs DF3
plot(siglda.res.climate$x[,2],
     siglda.res.climate$x[,3],
     ann=F,
     pch=sybs.climate,
     col=cols.climate,
     bg=bg.climate,
     las=1,
     frame=F,cex=2.5)

mtext(outer=F, side=1, line=2, text="DF2 (30.92%)",cex=1)
mtext(outer=F, side=2, line=2, text="DF3 (10.89%)",cex=1)
textxy(siglda.res.climate$x[,2], siglda.res.climate$x[,3], cex=1, labs = Community.proportion.climate$Community)

##Structure matrix plot CF1 and CF3
plot(S.climate[,c(2,3)], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:19){
  arrows(0,0,as.numeric(S.climate[i,2]),as.numeric(S.climate[i,3]), 
         lty=1, 
         lwd=1, 
         col="chocolate1",
         length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S.climate[,2]), as.numeric(S.climate[,3]),cex=1, labs=rownames(S.climate))
mtext(outer=F, side=1, line=2, text="DF2",las=1,cex=1)


#######################
## plots with environmental factors: precipitation, temperature, elevations
## use the ggplot2 to plot the data
gg.data.climate <- siglda.res.climate$x
gg.data.climate <- gg.data.climate[1:98,]
gg.data.climate <- cbind (gg.data.climate, Community.proportion.climate[,24:26])
gg.data.climate <- as.data.frame(gg.data.climate)

## DF1 against DF2
p1.climate <- ggplot(gg.data.climate, aes(LD1, LD2)) + 
  geom_point(aes(color=MAT), size = 3, alpha = 1) + 
  scale_colour_gradient(low ="#4169E1", high = "#EE5C42", breaks=c(24,16,8,0)) +
  theme(legend.position="none")
p2.climate<-ggplot(gg.data.climate, aes(LD1, LD2)) + 
  geom_point(aes(color=MAP), size = 3, alpha = 1) +
  scale_colour_gradient(low ="#E69F00" ,high = "#00FF7F", breaks=c(6000,4800,3600,2400,1200)) +
  theme(legend.position="none")
p3.climate<-ggplot(gg.data.climate, aes(LD1, LD2)) + 
  geom_point(aes(color=Elevation), size = 3, alpha = 1) +
  scale_color_gradient(low ="grey80", high = "grey10", breaks=c(3200, 2400, 1600, 800)) +
  theme(legend.position="none")
## DF1 against DF3
p4.climate <- ggplot(gg.data.climate, aes(LD1, LD3)) + 
  geom_point(aes(color=MAT), size = 3, alpha = 1) + 
  scale_colour_gradient(low ="#4169E1", high = "#EE5C42", breaks=c(24,16,8,0)) +
  theme(legend.position="none")
p5.climate<-ggplot(gg.data.climate, aes(LD1, LD3)) + 
  geom_point(aes(color=MAP), size = 3, alpha = 1) +
  scale_colour_gradient(low ="#E69F00" ,high = "#00FF7F", breaks=c(6000,4800,3600,2400,1200)) +
  theme(legend.position="none")
p6.climate<-ggplot(gg.data.climate, aes(LD1, LD3)) + 
  geom_point(aes(color=Elevation), size = 3, alpha = 1) +
  scale_color_gradient(low ="grey80", high = "grey10", breaks=c(3200, 2400, 1600, 800)) +
  theme(legend.position="none")
## DF2 against DF3
p7.climate <- ggplot(gg.data.climate, aes(LD2, LD3)) + 
  geom_point(aes(color=MAT), size = 3, alpha = 1) + 
  scale_colour_gradient(low ="#4169E1", high = "#EE5C42", breaks=c(24,16,8,0))
p8.climate<-ggplot(gg.data.climate, aes(LD2, LD3)) + 
  geom_point(aes(color=MAP), size = 3, alpha = 1) +
  scale_colour_gradient(low ="#E69F00" ,high = "#00FF7F", breaks=c(6000,4800,3600,2400,1200))
p9.climate<-ggplot(gg.data.climate, aes(LD2, LD3)) + 
  geom_point(aes(color=Elevation), size = 3, alpha = 1) +
  scale_color_gradient(low ="grey80", high = "grey10", breaks=c(3200, 2400, 1600, 800))

grid.arrange(p1.climate, p4.climate, p7.climate, 
             p2.climate, p5.climate, p8.climate, 
             p3.climate, p6.climate, p9.climate, 
             ncol=3, nrow=3)

###########################
## Eight Vegetation types
######################

## combine data
combined.data.vegetation<-rbind(Community.proportion.vegetation[,4:23], extinct.data[,4:23])

# number of each vegetation
n.tropical.rainforest<-length(which(Community.proportion.vegetation$Vegetation=="Tropical rainforest"))
n.tropical.seasonal.forest<-length(which(Community.proportion.vegetation$Vegetation=="Tropical seasonal forest"))
n.savanna<-length(which(Community.proportion.vegetation$Vegetation=="Savanna"))
n.grassland<-length(which(Community.proportion.vegetation$Vegetation=="Grassland"))
n.desert<-length(which(Community.proportion.vegetation$Vegetation=="Desert"))
n.shrubland<-length(which(Community.proportion.vegetation$Vegetation=="Shrubland"))
n.temperate.forest<-length(which(Community.proportion.vegetation$Vegetation=="Temperate forest"))
n.boreal.forest<-length(which(Community.proportion.vegetation$Vegetation=="Boreal forest"))

## Discriminant Functional Analysis on Significant Indices(<0.001)
sigdata.vegetation<-combined.data.vegetation

# Climate<-combined.data[,1]
Vegetation<-combined.data.vegetation[,1]

# ## data transformation transformed data = arcsin (sqrt(original data))
sigdata.vegetation<-cbind(Vegetation, asin(sqrt(sigdata.vegetation[,2:20])))

## apply X 100 to the precentage data
# sigdata<-cbind(Habitat, sigdata[,2:20]*100)

Strudata.vegetation<-sigdata.vegetation[1:n.community,2:20]

## Testing for significant differences between groups
sigindex.vegetation <- as.matrix(sigdata.vegetation[,c(2:20)])
sigoutput.vegetation <- manova(sigindex.vegetation~sigdata.vegetation$Vegetation)
sigsum.wilks.vegetation<-summary(sigoutput.vegetation, test="Wilks")
sigsum.avo.vegetation<-summary.aov(sigoutput.vegetation)

## DFA based vegetation types
siglda.mod.vegetation <- lda(Vegetation ~ ., 
                             data = sigdata.vegetation[1:n.community,], 
                             prior=c(1/8, 1/8, 1/8, 1/8, 1/8, 
                                     1/8, 1/8, 1/8))

siglda.res.vegetation <- predict(siglda.mod.vegetation, sigdata.vegetation)

## eigenvaules
siglda.mod.vegetation$svd

## First two coefficient of the LDA of the osteological indices
coemtr.vegetation<-cbind(siglda.mod.vegetation$scaling[,1],siglda.mod.vegetation$scaling[,2])

table_7vegetation<-table(sigdata.vegetation$Vegetation[1:n.community], siglda.res.vegetation$class[1:n.community])
write.table(table_7vegetation, file = "table_8vegetation.csv", sep = ",")

## grey color
cols.vegetation=c(rep("darkviolet", n.boreal.forest),         # boreal forest
                  rep("khaki", n.desert),         # desert
                  rep("goldenrod", n.grassland),        # grassland
                  rep("firebrick1", n.savanna),          # savanna
                  rep("darkgoldenrod1", n.shrubland),         # shrubland
                  rep("green3", n.temperate.forest),         # temperate forest 
                  rep("red", n.tropical.rainforest),         # tropical rainforest
                  rep("orangered", n.tropical.seasonal.forest),        # tropical seasonal forest
                  rep("grey30",5))    


# symbols based on the different habitats
sybs.vegetation=c(rep(1, n.boreal.forest),         # boreal forest
                  rep(18, n.desert),         # desert
                  rep(15, n.grassland),        # grassland
                  rep(2, n.savanna),          # savanna
                  rep(4, n.shrubland),         # shrubland
                  rep(0, n.temperate.forest),         # temperate forest 
                  rep(16, n.tropical.rainforest),         # tropical rainforest
                  rep(11, n.tropical.seasonal.forest),        # tropical seasonal forest
                  rep(16,5))    

## Structual Matrix of DFA 
StruCol.vegetation<-lda.structure(siglda.res.vegetation$x[1:n.community,],Strudata.vegetation) 
ld1.vegetation<-as.matrix(StruCol.vegetation$LD1)
ld2.vegetation<-as.matrix(StruCol.vegetation$LD2)
ld3.vegetation<-as.matrix(StruCol.vegetation$LD3)

write.csv(StruCol.vegetation, file = "StructMatrix_eight_vegetation.csv")
S.vegetation<-as.matrix(cbind(ld1.vegetation,ld2.vegetation,ld3.vegetation))

par(mfrow=c(3,2),oma=c(3,3,1,1),mar=c(2,2,0.5,1),
    mgp=c(3,0.5,0),tck=-0.005,cex=0.8,cex.axis=0.7)

## Plot DF1 vs DF2
plot(siglda.res.vegetation$x[,1],
     siglda.res.vegetation$x[,2],
     ann=F,
     pch=sybs.vegetation,
     col=cols.vegetation,
     las=1,
     frame=F,cex=2.5)

abline(v=0,lty=3)
abline(h=0,lty=3)
mtext(outer=F, side=2, line=2, text="DF2 (18.25%)",cex=1)
textxy(siglda.res.vegetation$x[,1], siglda.res.vegetation$x[,2], cex=1, 
       labs = Community.proportion.vegetation$Community)

## Structure matrix plot DF1 and DF2
plot(S.vegetation[,1:2], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:19) {
  arrows(0, 0, as.numeric(S.vegetation[i,1]), as.numeric(S.vegetation[i,2]), lty=1, lwd=1, 
         col="chocolate1", length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S.vegetation[,1]), as.numeric(S.vegetation[,2]), cex=1, labs=rownames(S.vegetation))
mtext(outer=F, side=3, line=2, text="Coefficient",las=1,cex=1)

## Plot DF1 vs DF3
plot(siglda.res.vegetation$x[,1],
     siglda.res.vegetation$x[,3],
     ann=F,
     pch=sybs.vegetation,
     col=cols.vegetation,
     las=1,
     frame=F,cex=2.5)

mtext(outer=F, side=1, line=1, text="DF1 (53.57%)",cex=1)
mtext(outer=F, side=2, line=2, text="DF3 (14.48%)",cex=1)
textxy(siglda.res.vegetation$x[,1], siglda.res.vegetation$x[,3], cex=1, 
       labs = Community.proportion.vegetation$Community)

##Structure matrix plot DF1 and DF3
plot(S.vegetation[,c(1,3)], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:19){
  arrows(0,0,as.numeric(S.vegetation[i,1]),as.numeric(S.vegetation[i,3]), 
         lty=1, 
         lwd=1, 
         col="chocolate1",
         length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S.vegetation[,1]), as.numeric(S.vegetation[,3]),cex=1, labs=rownames(S.vegetation))
mtext(outer=F, side=1, line=2, text="DF1",las=1,cex=1)


## Plot DF2 vs DF3
plot(siglda.res.vegetation$x[,2],
     siglda.res.vegetation$x[,3],
     ann=F,
     pch=sybs.vegetation,
     col=cols.vegetation,
     las=1,
     frame=F,cex=2.5)

mtext(outer=F, side=1.5, line=2, text="DF2 (18.25%)",cex=1)
mtext(outer=F, side=2, line=2, text="DF3 (14.48%)",cex=1)
textxy(siglda.res.vegetation$x[,2], siglda.res.vegetation$x[,3], cex=1, 
       labs = Community.proportion.vegetation$Community)


##Structure matrix plot DF2 and DF3
plot(S.vegetation[,c(2,3)], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:19){
  arrows(0,0,as.numeric(S.vegetation[i,2]),as.numeric(S.vegetation[i,3]), 
         lty=1, 
         lwd=1, 
         col="chocolate1",
         length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S.vegetation[,2]), as.numeric(S.vegetation[,3]),cex=1, labs=rownames(S.vegetation))
mtext(outer=F, side=1, line=2, text="DF2",las=1,cex=1)


### 3D plots
par(mfrow=c(1,1),oma=c(3,3,1,1),mar=c(0,1,0,1),
    mgp=c(3,0.5,0),tck=-0.005,cex=0.8,cex.axis=0.7)
s3d.vegetation<-scatterplot3d(siglda.res.vegetation$x[,1],
                   siglda.res.vegetation$x[,2],
                   siglda.res.vegetation$x[,3],
                   xlab="DF1(57.05%)",
                   ylab="DF2(17.94%)",
                   zlab="DF3(12.74%)",
                   pch=sybs.vegetation,color=cols.vegetation, angle = 60, type = "h")
s3d.coords.vegetation<-s3d.vegetation$xyz.convert(siglda.res.vegetation$x[,1],
                            siglda.res.vegetation$x[,2],
                            siglda.res.vegetation$x[,3])

text(x=s3d.coords.vegetation$x, 
     y=s3d.coords.vegetation$y,
     z=s3d.coords.vegetation$z,
     labels=Community.proportion.vegetation$Community,
     cex=0.5,pos=4, col = "red")

## 3D RGL plot
plot3d(siglda.res.vegetation$x[,1],
       siglda.res.vegetation$x[,2],
       siglda.res.vegetation$x[,3],
       pch=sybs.vegetation,
       col=cols.vegetation,
       size=10, 
       lwd=2,
       xlab = "",
       ylab = "",
       zlab = "")
decorate3d(axes = F,top = T)


######################################
## use the ggplot2 to plot the data
gg.data.vegetation <- siglda.res.vegetation$x
gg.data.vegetation <- gg.data.vegetation[1:98,]
gg.data.vegetation <- gg.data.vegetation[,1:3]
gg.data.vegetation <- cbind(gg.data.vegetation, Community.proportion.vegetation[,24:26])
gg.data.vegetation <- as.data.frame(gg.data.vegetation)

## DF1 against DF2
p1.vegetation <- ggplot(gg.data.vegetation, aes(LD1, LD2)) + 
  geom_point(aes(color=MAT), size = 3, alpha = 1) + 
  scale_colour_gradient(low ="#4169E1", high = "#EE5C42", breaks=c(24,16,8,0)) +
  theme(legend.position="none")
p2.vegetation<-ggplot(gg.data.vegetation, aes(LD1, LD2)) + 
  geom_point(aes(color=MAP), size = 3, alpha = 1) +
  scale_colour_gradient(low ="#E69F00" ,high = "#00FF7F", breaks=c(6000,4800,3600,2400,1200)) +
  theme(legend.position="none")
p3.vegetation<-ggplot(gg.data.vegetation, aes(LD1, LD2)) + 
  geom_point(aes(color=Elevation), size = 3, alpha = 1) +
  scale_color_gradient(low ="grey80", high = "grey10", breaks=c(3200, 2400, 1600, 800)) +
  theme(legend.position="none")
## DF1 against DF3
p4.vegetation <- ggplot(gg.data.vegetation, aes(LD1, LD3)) + 
  geom_point(aes(color=MAT), size = 3, alpha = 1) + 
  scale_colour_gradient(low ="#4169E1", high = "#EE5C42", breaks=c(24,16,8,0)) +
  theme(legend.position="none")
p5.vegetation<-ggplot(gg.data.vegetation, aes(LD1, LD3)) + 
  geom_point(aes(color=MAP), size = 3, alpha = 1) +
  scale_colour_gradient(low ="#E69F00" ,high = "#00FF7F", breaks=c(6000,4800,3600,2400,1200)) +
  theme(legend.position="none")
p6.vegetation<-ggplot(gg.data.vegetation, aes(LD1, LD3)) + 
  geom_point(aes(color=Elevation), size = 3, alpha = 1) +
  scale_color_gradient(low ="grey80", high = "grey10", breaks=c(3200, 2400, 1600, 800)) +
  theme(legend.position="none")
## DF2 against DF3
p7.vegetation <- ggplot(gg.data.vegetation, aes(LD2, LD3)) + 
  geom_point(aes(color=MAT), size = 3, alpha = 1) + 
  scale_colour_gradient(low ="#4169E1", high = "#EE5C42", breaks=c(24,16,8,0))
p8.vegetation<-ggplot(gg.data.vegetation, aes(LD2, LD3)) + 
  geom_point(aes(color=MAP), size = 3, alpha = 1) +
  scale_colour_gradient(low ="#E69F00" ,high = "#00FF7F", breaks=c(6000,4800,3600,2400,1200))
p9.vegetation<-ggplot(gg.data.vegetation, aes(LD2, LD3)) + 
  geom_point(aes(color=Elevation), size = 3, alpha = 1) +
  scale_color_gradient(low ="grey80", high = "grey10", breaks=c(3200, 2400, 1600, 800))

grid.arrange(p1.vegetation, p4.vegetation, p7.vegetation, 
             p2.vegetation, p5.vegetation, p8.vegetation, 
             p3.vegetation, p6.vegetation, p9.vegetation, 
             ncol=3, nrow=3)


###########################################################
## Regression against the DF1 of different environments
###########################################################
## Two habitats
## DF1 regression against environmental factors
## against MAP
DF1.habitat <- siglda.res.habitat$x[,1][1:98]
fit.habitat.precipitation <- lm(DF1.habitat~Community.proportion.habitat$MAP)
summary(fit.habitat.precipitation)
qqPlot(fit.habitat.precipitation)
leveragePlots(fit.habitat.precipitation)
## against Elevation
fit.habitat.elevation <- lm(DF1.habitat~Community.proportion.habitat$Elevation)
summary(fit.habitat.elevation)
qqPlot(fit.habitat.elevation)
leveragePlots(fit.habitat.elevation)
## against MAT
fit.habitat.temperature <- lm(DF1.habitat~Community.proportion.habitat$MAT)
summary(fit.habitat.temperature)
qqPlot(fit.habitat.temperature)
leveragePlots(fit.habitat.temperature)

## Four climates
## DF1 regression against environmental factors
## against MAP
DF1.climate <- siglda.res.climate$x[,1][1:98]
DF1.fit.climate.precipitation <- lm(DF1.climate~Community.proportion.climate$MAP)
summary(DF1.fit.climate.precipitation)
qqPlot(DF1.fit.climate.precipitation)
leveragePlots(DF1.fit.climate.precipitation)
## against Elevation
DF1.fit.climate.elevation <- lm(DF1.climate~Community.proportion.climate$Elevation)
summary(DF1.fit.climate.elevation)
qqPlot(DF1.fit.climate.elevation)
leveragePlots(DF1.fit.climate.elevation)
## against MAT
DF1.fit.climate.temperature <- lm(DF1.climate~Community.proportion.climate$MAT)
summary(DF1.fit.climate.temperature)
qqPlot(DF1.fit.climate.temperature)
leveragePlots(DF1.fit.climate.temperature)

## DF2 regression against environmental factors
## against MAP
DF2.climate <- siglda.res.climate$x[,2][1:98]
DF2.fit.climate.precipitation <- lm(DF2.climate~Community.proportion.climate$MAP)
summary(DF2.fit.climate.precipitation)
qqPlot(DF2.fit.climate.precipitation)
leveragePlots(DF2.fit.climate.precipitation)
## against Elevation
DF2.fit.climate.elevation <- lm(DF2.climate~Community.proportion.climate$Elevation)
summary(DF2.fit.climate.elevation)
qqPlot(DF2.fit.climate.elevation)
leveragePlots(DF2.fit.climate.elevation)
## against MAT
DF2.fit.climate.temperature <- lm(DF2.climate~Community.proportion.climate$MAT)
summary(DF2.fit.climate.temperature)
qqPlot(DF2.fit.climate.temperature)
leveragePlots(DF2.fit.climate.temperature)

## DF3 regression against environmental factors
## against MAP
DF3.climate <- siglda.res.climate$x[,3][1:98]
DF3.fit.climate.precipitation <- lm(DF3.climate~Community.proportion.climate$MAP)
summary(DF3.fit.climate.precipitation)
qqPlot(DF3.fit.climate.precipitation)
leveragePlots(DF2.fit.climate.precipitation)
## against Elevation
DF3.fit.climate.elevation <- lm(DF3.climate~Community.proportion.climate$Elevation)
summary(DF3.fit.climate.elevation)
qqPlot(DF3.fit.climate.elevation)
leveragePlots(DF3.fit.climate.elevation)
## against MAT
DF3.fit.climate.temperature <- lm(DF3.climate~Community.proportion.climate$MAT)
summary(DF3.fit.climate.temperature)
qqPlot(DF3.fit.climate.temperature)
leveragePlots(DF3.fit.climate.temperature)


## Eight vegetations
## DF1 regression against environmental factors
## against MAP
DF1.vegetation <- siglda.res.vegetation$x[,1][1:98]
fit.vegetation.precipitation <- lm(DF1.vegetation~Community.proportion.vegetation$MAP)
summary(fit.vegetation.precipitation)
qqPlot(fit.vegetation.precipitation)
leveragePlots(fit.vegetation.precipitation)
## against Elevation
fit.vegetation.elevation <- lm(DF1.vegetation~Community.proportion.vegetation$Elevation)
summary(fit.vegetation.elevation)
qqPlot(fit.vegetation.elevation)
leveragePlots(fit.vegetation.elevation)
## against MAT
fit.vegetation.temperature <- lm(DF1.vegetation~Community.proportion.vegetation$MAT)
summary(fit.vegetation.temperature)
qqPlot(fit.vegetation.temperature)
leveragePlots(fit.vegetation.temperature)

## DF2 regression against environmental factors
## against MAP
DF2.vegetation <- siglda.res.vegetation$x[,2][1:98]
DF2.fit.vegetation.precipitation <- lm(DF2.vegetation~Community.proportion.vegetation$MAP)
summary(DF2.fit.vegetation.precipitation)
qqPlot(DF2.fit.vegetation.precipitation)
leveragePlots(DF2.fit.vegetation.precipitation)
## against Elevation
DF2.fit.vegetation.elevation <- lm(DF2.vegetation~Community.proportion.vegetation$Elevation)
summary(DF2.fit.vegetation.elevation)
qqPlot(DF2.fit.vegetation.elevation)
leveragePlots(DF2.fit.vegetation.elevation)
## against MAT
DF2.fit.vegetation.temperature <- lm(DF2.vegetation~Community.proportion.vegetation$MAT)
summary(DF2.fit.vegetation.temperature)
qqPlot(DF2.fit.vegetation.temperature)
leveragePlots(DF2.fit.vegetation.temperature)

## DF3 regression against environmental factors
## against MAP
DF3.vegetation <- siglda.res.vegetation$x[,3][1:98]
DF3.fit.vegetation.precipitation <- lm(DF3.vegetation~Community.proportion.vegetation$MAP)
summary(DF3.fit.vegetation.precipitation)
qqPlot(DF3.fit.vegetation.precipitation)
leveragePlots(DF2.fit.vegetation.precipitation)
## against Elevation
DF3.fit.vegetation.elevation <- lm(DF3.vegetation~Community.proportion.vegetation$Elevation)
summary(DF3.fit.vegetation.elevation)
qqPlot(DF3.fit.vegetation.elevation)
leveragePlots(DF3.fit.vegetation.elevation)
## against MAT
DF3.fit.vegetation.temperature <- lm(DF3.vegetation~Community.proportion.vegetation$MAT)
summary(DF3.fit.vegetation.temperature)
qqPlot(DF3.fit.vegetation.temperature)
leveragePlots(DF3.fit.vegetation.temperature)

### all regression plots
par(mfrow=c(3,3),
    oma=c(3,3,1,1),
    mar=c(2,2,2.5,1),
    mgp=c(3,0.5,0),
    tck=-0.005,
    cex=0.8,
    cex.axis=0.7)

## two habitats
plot(Community.proportion.habitat$MAP, DF1.habitat, bty = "n")
abline(fit.habitat.precipitation, col="firebrick2", lwd=2)
r.squared.habitat.precipitation <- summary(fit.habitat.precipitation)$r.squared
habitat.precipitation.label <- bquote(italic(R)^2 == . (format(r.squared.habitat.precipitation, digits = 3)))
text(x=5200, y =-2.5, labels = habitat.precipitation.label)
mtext(outer=F, side=3, line=0, text="Precipitation against DF1", cex=1)
mtext(outer=F, side=2, line=2, text="Two habitats", cex=1)

plot(Community.proportion.habitat$Elevation, DF1.habitat, bty = "n")
abline(fit.habitat.elevation, col="firebrick2", lwd=2)
r.squared.habitat.elevation <- summary(fit.habitat.elevation)$r.squared
habitat.elevation.label <- bquote(italic(R)^2 == . (format(r.squared.habitat.elevation, digits = 3)))
text(x=3500, y =3, labels = habitat.elevation.label)
mtext(outer=F, side=3, line=0, text="Elevation against DF1", cex=1)

plot(Community.proportion.habitat$MAT, DF1.habitat, bty = "n")
abline(fit.habitat.temperature, col="firebrick2", lwd=2)
r.squared.habitat.temperature <- summary(fit.habitat.temperature)$r.squared
habitat.temperature.label <- bquote(italic(R)^2 == . (format(r.squared.habitat.temperature, digits = 3)))
text(x=-2, y =3, labels = habitat.temperature.label)
mtext(outer=F, side=3, line=0, text="Temperature against DF1", cex=1)


## four climates
plot(Community.proportion.climate$MAP, DF1.climate, bty = "n")
abline(DF1.fit.climate.precipitation, col = "firebrick2", lwd=2)
r.squared.climate.precipitation <- summary(DF1.fit.climate.precipitation)$r.squared
climate.precipitation.label <- bquote(italic(R)^2 == . (format(r.squared.climate.precipitation, digits = 3)))
text(x=5200, y =-3.5, labels = climate.precipitation.label)
mtext(outer=F, side=3, line=0, text="Precipitation against DF1", cex=1)
mtext(outer=F, side=2, line=2, text="Four climates", cex=1)

plot(Community.proportion.climate$Elevation, DF1.climate, bty = "n")
abline(fit.climate.elevation, col = "firebrick2", lwd=2)
r.squared.climate.elevation <- summary(fit.climate.elevation)$r.squared
climate.elevation.label <- bquote(italic(R)^2 == . (format(r.squared.climate.elevation, digits = 3)))
text(x=3500, y =4, labels = climate.elevation.label)
mtext(outer=F, side=3, line=0, text="Elevation against DF1", cex=1)

plot(Community.proportion.climate$MAT, DF1.climate, bty = "n")
abline(fit.climate.temperature, col="firebrick2", lwd=2)
r.squared.climate.temperature <- summary(fit.climate.temperature)$r.squared
climate.temperature.label <- bquote(italic(R)^2 == . (format(r.squared.climate.temperature, digits = 3)))
text(x=-2, y =4, labels = climate.temperature.label)
mtext(outer=F, side=3, line=0, text="Temperature against DF1", cex=1)


## eight vegetations
plot(Community.proportion.vegetation$MAP, DF1.vegetation, bty = "n")
abline(fit.vegetation.precipitation, col = "firebrick2", lwd=2)
r.squared.vegetation.precipitation <- summary(fit.vegetation.precipitation)$r.squared
vegetation.precipitation.label <- bquote(italic(R)^2 == . (format(r.squared.vegetation.precipitation, digits = 3)))
text(x=5200, y =-2, labels = vegetation.precipitation.label)
mtext(outer=F, side=3, line=0, text="Precipitation against DF1", cex=1)
mtext(outer=F, side=2, line=2, text="Eight vegetations", cex=1)

plot(Community.proportion.vegetation$Elevation, DF1.vegetation, bty = "n")
abline(fit.vegetation.elevation, col = "firebrick2", lwd=2)
r.squared.vegetation.elevation <- summary(fit.vegetation.elevation)$r.squared
vegetation.elevation.label <- bquote(italic(R)^2 == . (format(r.squared.vegetation.elevation, digits = 3)))
text(x=3500, y =5.5, labels = vegetation.elevation.label)
mtext(outer=F, side=3, line=0, text="Precipitation against DF1", cex=1)

plot(Community.proportion.vegetation$MAT, DF1.vegetation, bty = "n")
abline(fit.vegetation.temperature, col="firebrick2", lwd=2)
r.squared.vegetation.temperature <- summary(fit.vegetation.temperature)$r.squared
vegetation.temperature.label <- bquote(italic(R)^2 == . (format(r.squared.vegetation.temperature, digits = 3)))
text(x=-2, y =5.5, labels = vegetation.temperature.label)
mtext(outer=F, side=3, line=0, text="Temperature against DF1", cex=1)

## ggplots
## two habitats
habitat.data <- as.data.frame(cbind(Community.proportion.habitat$MAT, 
                                    Community.proportion.habitat$MAP, 
                                    Community.proportion.habitat$Elevation,
                                    DF1.habitat))
colnames(habitat.data) <-c("MAT", "MAP","Elevation","DF1")

## four climtes
climate.data <- as.data.frame(cbind(Community.proportion.climate$MAT,
                                    Community.proportion.climate$MAP, 
                                    Community.proportion.climate$Elevation,
                                    DF1.climate))
colnames(climate.data) <-c("MAT", "MAP","Elevation","DF1")

## eight vegetations
vegetation.data <- as.data.frame(cbind(Community.proportion.vegetation$MAT, 
                                       Community.proportion.vegetation$MAP,
                                       Community.proportion.vegetation$Elevation,
                                       DF1.vegetation))
colnames(vegetation.data) <-c("MAT", "MAP","Elevation","DF1")

## regression R and p values equation
lm_eqn <- function(y, x){
  m <- lm(y ~ x);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

## two habitats
ht <- ggplot(habitat.data, aes(x=MAT, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 20, y = 2.8, label = lm_eqn(habitat.data$DF1, habitat.data$MAT), parse = TRUE, size = 2.5)

hp <- ggplot(habitat.data, aes(x=MAP, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm)+
  geom_text(x = 4000, y = 2.8, label = lm_eqn(habitat.data$DF1, habitat.data$MAP), parse = TRUE, size = 2.5)

he <- ggplot(habitat.data, aes(x=Elevation, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 3000, y = -1.5, label = lm_eqn(habitat.data$DF1, habitat.data$Elevation), parse = TRUE, size = 2.5)

## four climates
ct <- ggplot(climate.data, aes(x=MAT, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 2, y = 3.75, label = lm_eqn(climate.data$DF1, climate.data$MAT), parse = TRUE, size = 2.5)

cp <- ggplot(climate.data, aes(x=MAP, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 1500, y = 6, label = lm_eqn(climate.data$DF1, climate.data$MAP), parse = TRUE, size = 2.5)

ce <- ggplot(climate.data, aes(x=Elevation, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 3000, y = 3, label = lm_eqn(climate.data$DF1, climate.data$Elevation), parse = TRUE, size = 2.5)

## eight vegetations
vt <- ggplot(vegetation.data, aes(x=MAT, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 2, y = 5, label = lm_eqn(vegetation.data$DF1, vegetation.data$MAT), parse = TRUE, size = 2.5)

vp <- ggplot(vegetation.data, aes(x=MAP, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 1800, y = 6.8, label = lm_eqn(vegetation.data$DF1, vegetation.data$MAP), parse = TRUE, size = 2.5)

ve <- ggplot(vegetation.data, aes(x=Elevation, y=DF1)) + 
  geom_point(shape=1) + 
  geom_smooth(method = lm) +
  geom_text(x = 3000, y = 4.5, label = lm_eqn(vegetation.data$DF1, vegetation.data$Elevation), parse = TRUE, size = 2.5)

### Run EDisp_ERich_DF1_Latitudeianl_plot codes first for regression.DF1.Latitude
grid.arrange(regression.DF1.latitude, ht, hp, he, ncol=2, nrow=2)

## 
grid.arrange(p1.climate, ct, p1.vegetation, vt, 
             p2.climate, cp, p2.vegetation, vp, 
             p3.climate, ce, p3.vegetation, ve, 
             ncol=4, nrow= 3)

#############################################################################################
### EDisp+ERich relationships with environemntal factors (latitude, MAT, MAP, and elevation)
### please run the EDisp_ERich_DF1_Latitudinal_plot R codes
### then use "coords.eco"
### environmental factor
#############################################################################################

setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis/Geographic Plots")

### GPS coordinates of each community
coords<-read.csv(file="GPS_coordinate_13July2017.csv", header=TRUE, sep = ",")

setwd("/Users/mengchen/Documents/Research Projects/4 Paleoecological studies/MM Community Project/Database and Analysis")

EDisp_ERich <- read.csv(file = "EDisp_ERich_30June2017.csv", sep = ",", header = T)

coords.eco.data <- cbind(coords, EDisp_ERich)
coords.eco.data <- as.data.frame(coords.eco.data)
## EDisp
Lat.EDisp <- ggplot(coords.eco.data, aes(abs(Lat), Mean.EDisp)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 50, y = 7, label = lm_eqn(coords.eco.data$Mean.EDisp, abs(coords.eco.data$Lat)), parse = TRUE, size = 3)
MAT.EDisp <- ggplot(coords.eco.data, aes(MAT, Mean.EDisp)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 2, y = 7, label = lm_eqn(coords.eco.data$Mean.EDisp, coords.eco.data$MAT), parse = TRUE, size = 3)
MAP.EDisp <- ggplot(coords.eco.data, aes(MAP, Mean.EDisp)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 4500, y = 7, label = lm_eqn(coords.eco.data$Mean.EDisp, coords.eco.data$MAP), parse = TRUE, size = 3)
Elevation.EDisp <- ggplot(coords.eco.data, aes(Elevation, Mean.EDisp)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 3000, y = 7, label = lm_eqn(coords.eco.data$Mean.EDisp, coords.eco.data$Elevation), parse = TRUE, size = 3)

## ERich
Lat.ERich <- ggplot(coords.eco.data, aes(abs(Lat), ERich)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 50, y = 20, label = lm_eqn(coords.eco.data$ERich, abs(coords.eco.data$Lat)), parse = TRUE, size = 3)
MAT.ERich <- ggplot(coords.eco.data, aes(MAT, ERich)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 2, y = 20, label = lm_eqn(coords.eco.data$ERich, coords.eco.data$MAT), parse = TRUE, size = 3)
MAP.ERich <- ggplot(coords.eco.data, aes(MAP, ERich)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 4500, y = 20, label = lm_eqn(coords.eco.data$ERich, coords.eco.data$MAP), parse = TRUE, size = 3)
Elevation.ERich <- ggplot(coords.eco.data, aes(Elevation, ERich)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 3000, y = 20, label = lm_eqn(coords.eco.data$ERich, coords.eco.data$Elevation), parse = TRUE, size = 3)

## MAT
ef.factor.MAT <- ggplot(coords.eco.data, aes(abs(Lat), MAT)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 50, y = 30, label = lm_eqn(coords.eco.data$MAT, abs(coords.eco.data$Lat)), parse = TRUE, size = 2.5)
### MAP
ef.factor.MAP <- ggplot(coords.eco.data, aes(abs(Lat), MAP)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 50, y = 5000, label = lm_eqn(coords.eco.data$MAP, abs(coords.eco.data$Lat)), parse = TRUE, size = 2.5)
### Elevation
ef.factor.Elevation <- ggplot(coords.eco.data, aes(abs(Lat), Elevation)) +
  geom_point(shape=1) +
  geom_smooth(method = lm) +
  geom_text(x = 50, y = 3000, label = lm_eqn(coords.eco.data$Elevation, abs(coords.eco.data$Lat)), parse = TRUE, size = 2.5)


grid.arrange(ef.factor.MAT, p1.climate, ct, p1.vegetation, vt, 
             ef.factor.MAP, p2.climate, cp, p2.vegetation, vp, 
             ef.factor.Elevation, p3.climate, ce, p3.vegetation, ve, 
             ncol=5, nrow= 3)

grid.arrange(Lat.EDisp, Lat.ERich, ncol = 2)

grid.arrange(Lat.EDisp, MAT.EDisp, MAP.EDisp, Elevation.EDisp,
             Lat.ERich, MAT.ERich, MAP.ERich, Elevation.ERich,
             ncol = 4, nrow=2)

#################################################################
### predict the paleclimates of exticnt mammalian communities
#################################################################

## functions based on paleolatitude
MAT.lm.latitude <- lm(MAT~abs(Lat), data = coords.eco.data)
MAP.lm.latitude <- lm(MAP~abs(Lat), data = coords.eco.data)
## two habitats function
MAT.lm.habitat <- lm(MAT~DF1.habitat, data = habitat.data)
MAP.lm.habitat <- lm(MAP~DF1.habitat, data = habitat.data)
# Elevation.lm.habitat <- lm(Elevation~DF1.habitat, data = habitat.data)
## four climates function
MAT.lm.climate <- lm(MAT~DF1.climate, data = climate.data)
MAP.lm.climate <- lm(MAP~DF1.climate, data = climate.data)
# Elevation.lm.climate <- lm(Elevation~DF1.climate, data = climate.data)
## eight climates function
MAT.lm.vegetation <- lm(MAT~DF1.vegetation, data = vegetation.data)
MAP.lm.vegetation <- lm(MAP~DF1.vegetation, data = vegetation.data)
# Elevation.lm.vegetation <- lm(Elevation~DF1.vegetation, data = vegetation.data)

## Prediciton based on the latitude
## Cretaceous
paleolatitude.K <- c(41.9, 35.3, 48.5)  # of Mesozic mammalian communtiies
paleolatitude.K <- as.data.frame(paleolatitude.K)
rownames <- c("mean", "lwr",  "upr")
colnames(paleolatitude.K) <- c("Lat")
## prdiction
MAT.prediction.latitude.K <- predict(MAT.lm.latitude, paleolatitude.K, interval = "confidence", 
                                   level = 0.95)
rownames(MAT.prediction.latitude.K) <- c("Mean", "Lower", "Upper")
colnames(MAT.prediction.latitude.K) <- c("fit.MAT", "lwr.MAT", "upr.MAT")

MAP.prediction.latitude.K <- predict(MAP.lm.latitude, paleolatitude.K, interval = "confidence", 
                                   level = 0.95)
rownames(MAP.prediction.latitude.K) <- c("Mean", "Lower", "Upper")
colnames(MAP.prediction.latitude.K) <- c("fit.MAP", "lwr.MAP", "upr.MAP")

prediction.latitude.K <- cbind(MAT.prediction.latitude.K, MAP.prediction.latitude.K)
prediction.latitude.K <- as.data.frame(prediction.latitude.K)
prediction.latitude.K$community <- c("K.Mean")

## data of DF1 of two habitat DFA
DF1.extinct.habitat <- siglda.res.habitat$x[99:103,]
## data of DF1 of four climate DFA
DF1.extinct.climate <- siglda.res.climate$x[99:103,][,1]
## data of DF1 of eight vegetation DFA
DF1.extinct.vegetation <- siglda.res.vegetation$x[99:103,][,1]

## exticnt DF1
DF1.extinct.all <- rbind(DF1.extinct.habitat, DF1.extinct.climate, DF1.extinct.vegetation)
DF1.extinct.all <- t (DF1.extinct.all)
rownames(DF1.extinct.all) <- extinct.data$Commuity
colnames(DF1.extinct.all) <- c("DF1.habitat", "DF1.climate", "DF1.vegetation")

MAT.prediction.habitat <- data.frame()
MAP.prediction.habitat <- data.frame()
# Elevation.prediction.habitat <- data.frame()

### for 2 habitats
for (i in 1:5) {
    new.data <- data.frame(DF1.habitat = DF1.extinct.all[i,1])
    
    temp.MAT.prediction <- predict(MAT.lm.habitat, 
                                   as.data.frame(new.data), 
                                   interval = "confidence", 
                                   level = 0.95)
    MAT.prediction.habitat <- rbind(MAT.prediction.habitat, c(temp.MAT.prediction))
    
    temp.MAP.prediction <- predict(MAP.lm.habitat, 
                                   as.data.frame(new.data), 
                                   interval = "confidence", 
                                   level = 0.95)
    MAP.prediction.habitat <- rbind(MAP.prediction.habitat, c(temp.MAP.prediction))
    
    # temp.Elevation.prediction <- predict(Elevation.lm.habitat, 
    #                                      as.data.frame(new.data), 
    #                                      interval = "confidence", 
    #                                      level = 0.95)
    # Elevation.prediction.habitat <- rbind(Elevation.prediction.habitat, c(temp.Elevation.prediction))
}

colnames(MAT.prediction.habitat) <- c("fit.MAT", "lwr.MAT", "upr.MAT")
rownames(MAT.prediction.habitat) <- extinct.data$Commuity
colnames(MAP.prediction.habitat) <- c("fit.MAP", "lwr.MAP", "upr.MAP")
rownames(MAP.prediction.habitat) <- extinct.data$Commuity

prediction.habitat <- cbind(MAT.prediction.habitat, MAP.prediction.habitat)
prediction.habitat$community <- extinct.data$Commuity
prediction.habitat <- rbind(prediction.habitat, prediction.latitude.K[1,])
prediction.habitat <- as.data.frame(prediction.habitat)
write.csv(prediction.habitat, file = "Results_Predicted_MAP_MAT_TwoHabitats.csv")

### for four climates
MAT.prediction.climate <- data.frame()
MAP.prediction.climate <- data.frame()

for (i in 1:5) {
  new.data <- data.frame(DF1.climate = DF1.extinct.all[i,2])
  temp.MAT.prediction <- predict(MAT.lm.climate, 
                                 as.data.frame(new.data), 
                                 interval = "confidence", 
                                 level = 0.95)
  MAT.prediction.climate <- rbind(MAT.prediction.climate, c(temp.MAT.prediction))
  
  temp.MAP.prediction <- predict(MAP.lm.climate, 
                                 as.data.frame(new.data), 
                                 interval = "confidence", 
                                 level = 0.95)
  MAP.prediction.climate<- rbind(MAP.prediction.climate, c(temp.MAP.prediction))
  
  # temp.Elevation.prediction <- predict(Elevation.lm.habitat, 
  #                                      as.data.frame(new.data), 
  #                                      interval = "confidence", 
  #                                      level = 0.95)
  # Elevation.prediction.habitat <- rbind(Elevation.prediction.habitat, c(temp.Elevation.prediction))
}

colnames(MAT.prediction.climate) <- c("fit.MAT", "lwr.MAT", "upr.MAT")
rownames(MAT.prediction.climate) <- extinct.data$Commuity
colnames(MAP.prediction.climate) <- c("fit.MAP", "lwr.MAP", "upr.MAP")
rownames(MAP.prediction.climate) <- extinct.data$Commuity

prediction.climate <- cbind(MAT.prediction.climate, MAP.prediction.climate)
prediction.climate$community <- extinct.data$Commuity
prediction.climate <- rbind(prediction.climate, prediction.latitude.K[1,])
prediction.climate <- as.data.frame(prediction.climate)
write.csv(prediction.climate, file = "Results_Predicted_MAP_MAT_FourClimates.csv")

### eight vegetations
MAT.prediction.vegetation <- data.frame()
MAP.prediction.vegetation <- data.frame()

for (i in 1:5) {
  new.data <- data.frame(DF1.vegetation = DF1.extinct.all[i,2])
  temp.MAT.prediction <- predict(MAT.lm.vegetation, 
                                 as.data.frame(new.data), 
                                 interval = "confidence", 
                                 level = 0.95)
  MAT.prediction.vegetation <- rbind(MAT.prediction.vegetation, c(temp.MAT.prediction))
  
  temp.MAP.prediction <- predict(MAP.lm.vegetation, 
                                 as.data.frame(new.data), 
                                 interval = "confidence", 
                                 level = 0.95)
  MAP.prediction.vegetation<- rbind(MAP.prediction.vegetation, c(temp.MAP.prediction))
  
  # temp.Elevation.prediction <- predict(Elevation.lm.habitat, 
  #                                      as.data.frame(new.data), 
  #                                      interval = "confidence", 
  #                                      level = 0.95)
  # Elevation.prediction.habitat <- rbind(Elevation.prediction.habitat, c(temp.Elevation.prediction))
}

colnames(MAT.prediction.vegetation) <- c("fit.MAT", "lwr.MAT", "upr.MAT")
rownames(MAT.prediction.vegetation) <- extinct.data$Commuity
colnames(MAP.prediction.vegetation) <- c("fit.MAP", "lwr.MAP", "upr.MAP")
rownames(MAP.prediction.vegetation) <- extinct.data$Commuity

prediction.vegetation <- cbind(MAT.prediction.vegetation, MAP.prediction.vegetation)
prediction.vegetation$community <- extinct.data$Commuity
prediction.vegetation <- rbind(prediction.vegetation, prediction.latitude.K[1,])
prediction.vegetation <- as.data.frame(prediction.vegetation)
write.csv(prediction.vegetation, file = "Results_Predicted_MAP_MAT_EightVegetations.csv")

## Whittacker plots of extinct mammalian communities
WhittakerPlot_Extinct_Habitat <- ggplot(data=prediction.habitat[-6,], aes(x=fit.MAT, y=fit.MAP, group=community)) + 
  geom_errorbar(aes(ymin=lwr.MAP, ymax=upr.MAP, color=community), width = 2) +
  geom_errorbarh(aes(xmin=lwr.MAT, xmax=upr.MAT, color=community)) +
  geom_point(aes(shape=community, color=community), size =4) +
  scale_shape_manual(values = c(15, 16, 17, 9, 18)) +
  scale_color_manual(values = c("#f8766d","#d7d78f", "#4fd3a5", "#00b0f6", "#ed90f6")) +
  coord_cartesian(xlim = c(-10, 35), ylim = c(0, 4500)) +
  scale_x_reverse() +
  theme(legend.position="none")

WhittakerPlot_Extinct_Climate<- ggplot(data=prediction.climate[-6,], aes(x=fit.MAT, y=fit.MAP, group=community)) + 
  geom_errorbar(aes(ymin=lwr.MAP, ymax=upr.MAP, color=community), width = 2) +
  geom_errorbarh(aes(xmin=lwr.MAT, xmax=upr.MAT, color=community)) +
  geom_point(aes(shape=community, color=community), size =4) +
  scale_shape_manual(values = c(15, 16, 17, 9, 18)) +
  scale_color_manual(values = c("#f8766d","#d7d78f", "#4fd3a5", "#00b0f6", "#ed90f6")) +
  coord_cartesian(xlim = c(-10, 35), ylim = c(0, 4500)) +
  scale_x_reverse() +
  theme(legend.position="none")

WhittakerPlot_Extinct_Vegetation <- ggplot(data=prediction.vegetation[-6,], aes(x=fit.MAT, y=fit.MAP, group=community)) + 
  geom_errorbar(aes(ymin=lwr.MAP, ymax=upr.MAP, color=community), width = 2) +
  geom_errorbarh(aes(xmin=lwr.MAT, xmax=upr.MAT, color=community)) +
  geom_point(aes(shape=community, color=community), size =4) +
  scale_shape_manual(values = c(15, 16, 17, 9, 18)) +
  scale_color_manual(values = c("#f8766d","#d7d78f", "#4fd3a5", "#00b0f6", "#ed90f6")) +
  coord_cartesian(xlim = c(-10, 35), ylim = c(0, 4500)) +
  scale_x_reverse()  +
  theme(legend.position="none")


## Whittaker plot of extant small-bodied mammlaina commmunities
## based on two haobitats
WhittakerPlot_habitat <- ggplot(data=coords.eco.data, aes(x=MAT, y=MAP, group=Habitat)) + 
  geom_point(aes(shape=Habitat, color=Habitat), size =4) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = c("grey30", "black")) +
  coord_cartesian(xlim = c(-10, 35), ylim = c(0, 4500)) +
  scale_x_reverse()+
  theme(legend.position="none")

## based on four climates
## the shape of the points have been changed to solid points in corresponding to the original ones
WhittakerPlot_climate <- ggplot(data=coords.eco.data, aes(x=MAT, y=MAP, group=Climate)) + 
  geom_point(aes(shape=Climate, color=Climate), size =4) +
  scale_shape_manual(values = c(18, 15, 17, 16)) +
  scale_color_manual(values = c("#59505A","dodgerblue3", "green3", "darkorange")) +
  coord_cartesian(xlim = c(-10, 35), ylim = c(0, 4500)) +
  scale_x_reverse()+
  theme(legend.position="none")

## based on eight vegetations
WhittakerPlot_vegetation<- ggplot(data=coords.eco.data, aes(x=MAT, y=MAP, group=Vegetation)) + 
  geom_point(aes(shape=Vegetation, color=Vegetation), size =4) +
  scale_shape_manual(values = c(1, 18, 15, 2, 4, 0, 16, 11)) +
  scale_color_manual(values = c("darkviolet","khaki", "goldenrod", "firebrick1","darkgoldenrod1", "green3", "red", "orangered")) +
  coord_cartesian(xlim = c(-10, 35), ylim = c(0, 4500)) +
  scale_x_reverse()+
  theme(legend.position="none")


grid.arrange(WhittakerPlot_habitat, WhittakerPlot_climate, WhittakerPlot_vegetation,
             WhittakerPlot_Extinct_Habitat, WhittakerPlot_Extinct_Climate, WhittakerPlot_Extinct_Vegetation,
             ncol=3, nrow= 2)
