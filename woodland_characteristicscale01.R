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

#set study area
studyExtent<-c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min x, max x, min y, max y
#now crop our points to this area
C<-crop(melesmeles.sp,studyExtent)
#read raster data
LCM=rast("LCMUK.tif")
melesmelesFin<-project(C,crs(LCM))
melesmelesCoords<-crds(melesmelesFin)

#Select a larger window
x.min <- min(melesmelesCoords[,1]) - 5000
x.max <- max(melesmelesCoords[,1]) + 5000
y.min <- min(melesmelesCoords[,2]) - 5000
y.max <- max(melesmelesCoords[,2]) + 5000
extent.new <- ext(x.min, x.max, y.min, y.max)
LCM <- crop(LCM$LCMUK_1, extent.new)


plot(LCM)
plot(melesmelesFin,add=TRUE)

#create a bakcground spatialPoints layer from the back.xy matrix
set.seed(11)
back.xy <- spatSample(LCM, size=1000,as.points=TRUE) 
plot(LCM)
plot(melesmelesFin,add=T)
plot(back.xy,add=TRUE, col='red')

#Extract values from the background points
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

#create a data frame of absence/background point coordinates
Abs<-data.frame(crds(back.xy),Pres=0)

#inspect the first few rows of the absence/background dataset
head(Abs)

Pres<-data.frame(crds(melesmelesFin),Pres=1)

#inspect the first few rows of the presence dataset
head(Pres)


#bind the two data frames by row 
melesmelesData<-rbind(Pres,Abs)


melesmelesSF=st_as_sf(melesmelesData,coords=c("x","y"),crs="EPSG:27700")

#access levels of the raster by treating them as categorical data ('factors' in R)
LCM<-as.factor(LCM)
levels(LCM)

#create an vector object called reclass
reclass <- c(0,1,rep(0,19))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix<- cbind(levels(LCM)[[1]],reclass)
RCmatrix<-RCmatrix[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)
#Use the reclassify() function to asssign new values to LCM with our reclassification matrix
RCmatrix

broadleaf <- classify(LCM, RCmatrix)

#reset the plotting panel
dev.off()

#plot
plot(broadleaf)
plot(melesmelesFin,add=TRUE)