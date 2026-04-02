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
#reset the plotting panel
dev.off()


#create an object to hold the distance parameter for  buffer
buf5km<-5000 

#use the st_buffer() function from the sf package applied to the first item of melesmelesFin

buffer.site1.5km<-st_buffer(melesmelesSF[1,],dist=buf5km) 
zoom(urban,buffer.site1.5km) # use the zoom() function for a close-up of the result.

plot(buffer.site1.5km$geometry,border="red",lwd=2,add=T) # add the buffer

#crop the urban layer to the extent of the buffer
buffer5km <- crop(urban, buffer.site1.5km) 
#clip the above again to the circle described by the buffer (doing this speeds up the process compared to using only the mask() function)
bufferlandcover5km <- mask(urban, buffer.site1.5km)

#calculate the area of the buffer according to the buffer width
bufferArea <- (3.14159*buf5km^2) 

#total area of woodland inside the buffer (625 is the area in metres of each cell of our 25x25m raster) na.rm=T makes sure any NA values are removed from the calculation (otherwise NA is returned) 
landcover5km <- sum(values(bufferlandcover5km),na.rm=T)*625 

#calculate percentage
percentlandcover5km <- landcover5km/bufferArea*100 

#return the result
percentlandcover5km



#function for automating whole dataset. The function is set up take two arguments: a data frame and a series of buffer distances.
landBuffer <- function(speciesData, r){         
  
  #buffer each point
  melesmelesBuffer <- st_buffer(speciesData, dist=r)                     
  
  #crop the woodland layer to the buffer extent
  bufferlandcover <- crop(urban, melesmelesBuffer)              
  
  # now extract the raster values (which should all be 1 for woodland and 0 for everything else) within each buffer and sum to get number of woodland cells inside the buffers.
  masklandcover <- extract(bufferlandcover, melesmelesBuffer,fun="sum")      
  #get woodland area (625 is the area in metres of each cell of our 25m raster)
  landcoverArea <- masklandcover$LCMUK_1*625  
  
  # convert to precentage cover (we use the st_area() function from the sf package to get the area of our buffer) but convert to a numeric object (because sf applies units i.e. metres which then cant be entered into numeric calculations)
  percentcover <- landcoverArea/as.numeric(st_area(melesmelesBuffer))*100 
  
  # return the result
  return(percentcover)                                       
}



