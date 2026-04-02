setwd("D:/Project/SE/SpatialEcology_Assessment1")

#Install packages
install.packages(c("terra","sf","mapview"))
library(terra)
library(sf)


#read data
melesmeles<- read.csv("Melesmeles.csv")
#data cleaning to ensure there are no NA
#subset the data to only include points with complete coordinates
melesmeles<-melesmeles[!is.na(melesmeles$Latitude),]
#remove all points with uncertainty > 1000m
melesmeles<-melesmeles[melesmeles$Coordinate.uncertainty_m<1001,]


#make spatial points layer
#create crs object
melesmeles.latlong<-data.frame(x=melesmeles$Longitude,y=melesmeles$Latitude)
#Use coordinates object to create our spatial points object
melesmeles.sp<-vect(melesmeles.latlong,geom=c("x","y"))
#check that the points now have our desired crs. 
crs(melesmeles.sp)<-"epsg:4326"

#show map
plot(melesmeles.sp)
