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


#Generate pseudo-absence points
set.seed(11)
back.xy <- spatSample(LCM, size=1000,as.points=TRUE) 
#create a spatialPoints layer from the back.xy matrix
plot(melesmelesFin,add=T)
plot(back.xy,add=TRUE, col='red')

#extract land cover values at background and presence points
eA<-extract(LCM,back.xy)
eP<-extract(LCM,melesmelesFin)

#calculate point frequency per landcover category (headings are landcover category codes - see LCM documentation for descriptions)
table(eA[,2])
table(eP[,2])

par(mfrow=c(1,2))

hist(eA[,2],freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,1))
hist(eP[,2],freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,1))

hist(eA[,2],freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,0.4))
hist(eP[,2],freq=FALSE,breaks=c(0:21),xlim=c(0,21),ylim=c(0,0.4))


#Build presence/absence dataset
Abs<-data.frame(crds(back.xy),Pres=0)
head(Abs)

Pres<-data.frame(crds(melesmelesFin),Pres=1)
head(Pres)

#bind the two data frames by row (both dataframes have the same column headings)
melesmelesData<-rbind(Pres,Abs)
#inspect
head(melesmelesData)

#convert to sf object for buffering
melesmelesSF=st_as_sf(melesmelesData,coords=c("x","y"),crs="EPSG:27700")

#access levels of the raster by treating them as categorical data ('factors' in R)
#Reclassify

#treat raster as categorical data
LCM <- as.factor(LCM)

#default = 0 for all classes
reclass <- rep(0, nrow(levels(LCM)[[1]]))

#choose urban classes
urban_codes <- c(20, 21)

#set urban/suburban to 1
reclass[levels(LCM)[[1]]$ID %in% urban_codes] <- 1

# create reclassification matrix
RCmatrix <- cbind(levels(LCM)[[1]]$ID, reclass)
RCmatrix <- apply(RCmatrix, 2, as.numeric)

# classify raster
urban <- classify(LCM, RCmatrix)

plot(urban)
plot(melesmelesFin, add = TRUE)




