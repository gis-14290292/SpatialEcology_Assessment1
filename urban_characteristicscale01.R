setwd("D:/Project/SE/SpatialEcology_Assessment1")

#Install packages
install.packages(c("terra","sf","mapview"))
library(terra)
library(sf)

#read data and clean
melesmeles<- read.csv("Melesmeles.csv")
head(melesmeles) 
class(melesmeles)

#data cleaning to ensure there are no NA
melesmeles<-melesmeles[!is.na(melesmeles$Latitude),]
#remove all points with uncertainty > 1000m
melesmeles<-melesmeles[melesmeles$Coordinate.uncertainty_m<1001,]
#create crs object
melesmeles.latlong<-data.frame(x=melesmeles$Longitude,
                               y=melesmeles$Latitude)
#Use coordinates object to create our spatial points object
melesmeles.sp<-vect(melesmeles.latlong,geom=c("x","y"))
#check that the points now have our desired crs. 
crs(melesmeles.sp)<-"epsg:4326"

plot(melesmeles.sp)
#define study area extent
studyExtent<-c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min x, max x, min y, max y

#now crop points to this area
C<-crop(melesmeles.sp,studyExtent)
#read in raster data
LCM=rast("LCMUK.tif")

#project points to same CRS as raster
melesmelesFin<-project(C,crs(LCM))
melesmelesCoords<-crds(melesmelesFin)

#make a slightly larger extent around points
x.min <- min(melesmelesCoords[,1]) - 5000
x.max <- max(melesmelesCoords[,1]) + 5000
y.min <- min(melesmelesCoords[,2]) - 5000
y.max <- max(melesmelesCoords[,2]) + 5000
extent.new <- ext(x.min, x.max, y.min, y.max)

# crop raster to this extent
LCM <- crop(LCM$LCMUK_1, extent.new)

plot(LCM)
plot(melesmelesFin,add=TRUE)