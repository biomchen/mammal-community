#######################################################
#### This function is used for calculate the EDisp ####
#### of speices within same commmunities ##############
#######################################################

EDisp_func_habitat <- function(data.input, habitat) {
    data.eds.mean<-data.frame()
    if (habitat == "Open") {
      data.input<-data.input[data.input$Habitat == "Open",]
    } 
    if (habitat == "Close") {
      data.input<-data.input[data.input$Habitat == "Close",]
    } 
    ### communities select for certain environment
    m<-unique(data.input$Community)
    for (n in m) {
      data.temp<-data.input[data.input$Community==n,]
      data.select.community<-data.temp[c(3,4,5)] ### the function trait data of the community selected
      for (i in 1:(dim(data.select.community)[1]-1)) { ### select first species to penultimate species within the community
        for (j in i:(dim(data.select.community)[1]-1)) {### select second to the last species with the communitiy
          data.eds.abs<-abs(data.select.community[i,]-data.select.community[j+1,]) ## absoluate different betweeen each function trait
          data.eds.pair<-rowSums(data.eds.abs) ### ecological disparity between a pair of species
          data.eds<-rbind(data.eds,c(data.eds.pair)) ### adding all results of all communities 
         }
      }
      data.eds.mean<-rbind(data.eds.mean, c(mean(data.eds[,1])))
      }
  return (data.eds.mean)
}

EDisp_func_climate <- function(data.input, climate) {
  data.eds.mean<-data.frame()
  data.eds<-data.frame()
  if (climate == "Tropical") {
    data.input<-data.input[data.input$Climate == "Tropical",]
  } 
  if (climate == "Arid") {
    data.input<-data.input[data.input$Climate == "Arid",]
  }
  if (climate == "Temperate") {
    data.input<-data.input[data.input$Climate == "Temperate",]
  } 
  if (climate == "Cold") {
    data.input<-data.input[data.input$Climate == "Cold",]
  } 
  ### communities select for certain environment
  m<-unique(data.input$Community)
  for (n in m) {
    data.temp<-data.input[data.input$Community==n,]
    data.select.community<-data.temp[c(3,4,5)] ### the function trait data of the community selected
    for (i in 1:(dim(data.select.community)[1]-1)) { ### select first species to penultimate species within the community
      for (j in i:(dim(data.select.community)[1]-1)) {### select second to the last species with the communitiy
        data.eds.abs<-abs(data.select.community[i,]-data.select.community[j+1,]) ## absoluate different betweeen each function trait
        data.eds.pair<-rowSums(data.eds.abs) ### ecological disparity between a pair of species
        data.eds<-rbind(data.eds,c(data.eds.pair)) ### adding all results of all communities 
      }
    }
    data.eds.mean<-rbind(data.eds.mean, c(mean(data.eds[,1])))
  }
  return (data.eds.mean)
}

EDisp_func_vegetation <- function(data.input, vegetation) {
  data.eds.mean<-data.frame()
  if (vegetation == "Tropical rainforest") {
    data.input<-data.input[data.input$Vegetation == "Tropical rainforest",]
  } 
  if (vegetation == "Tropical seasonal forest") {
    data.input<-data.input[data.input$Vegetation == "Tropical seasonal forest",]
  }
  if (vegetation == "Temperate forest") {
    data.input<-data.input[data.input$Vegetation == "Temperate forest",]
  }
  if (vegetation == "Grassland") {
    data.input<-data.input[data.input$Vegetation == "Grassland",]
  }
  if (vegetation == "Desert") {
    data.input<-data.input[data.input$Vegetation == "Temperate desert",]
  } 
  if (vegetation == "Savanna") {
    data.input<-data.input[data.input$Vegetation == "Savanna",]
  }
  if (vegetation == "Shrubland") {
    data.input<-data.input[data.input$Vegetation == "Shrubland",]
  } 
  if (vegetation == "Boreal forest") {
    data.input<-data.input[data.input$Vegetation == "Boreal forest",]
  } 
  ### communities select for certain environment
  m<-unique(data.input$Community)
  for (n in m) {
    data.temp<-data.input[data.input$Community==n,]
    data.select.community<-data.temp[c(3,4,5)] ### the function trait data of the community selected
    for (i in 1:(dim(data.select.community)[1]-1)) { ### select first species to penultimate species within the community
      for (j in i:(dim(data.select.community)[1]-1)) {### select second to the last species with the communitiy
        data.eds.abs<-abs(data.select.community[i,]-data.select.community[j+1,]) ## absoluate different betweeen each function trait
        data.eds.pair<-rowSums(data.eds.abs) ### ecological disparity between a pair of species
        data.eds<-rbind(data.eds,c(data.eds.pair)) ### adding all results of all communities 
      }
    }
    data.eds.mean<-rbind(data.eds.mean, c(mean(data.eds[,1])))
  }
  return (data.eds.mean)
}