#Install packages
install.packages(c("dismo","terra","spatstat","cowplot","ggplot2","precrec","glmnet","maxnet","ranger"))
#We also need to load the "MLR" package for this practical
install.packages("mlr")
#Function to convert raster to images for spatstat. Takes one argument "im" that should be a raster object
install.packages(c("terra","sf","mapview"))
library(terra)
library(sf)
library(spatstat) #for point process modelling and converting between raster and pixel image objectslibrary(spatstat)


#for point process modelling and converting between raster and pixel image objects
raster.as.im = function(im) {
  #get the resolution (cell size of the raster)
  r = raster::res(im)[1]
  #get the origin (bottom left corner of the raster/image)
  orig = ext(im)[c(1,3)]
  #set the coordinates of the columns which is just a series of number from zero increasing by 100 metres (the resolution of the raster) for every cell along the rows and columns.
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r)
  #set the coordinates of the columns
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r)
  
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ]
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}



#Read and clean species data
melesmeles<- read.csv("Melesmeles.csv")

#data cleaning to ensure there are no NA
melesmeles<-melesmeles[!is.na(melesmeles$Latitude),]

#remove all points with uncertainty > 1000m
melesmeles<-melesmeles[melesmeles$Coordinate.uncertainty_m<1001,]

#create crs object
melesmeles.latlong<-data.frame(x=melesmeles$Longitude,y=melesmeles$Latitude)

#Use coordinates object to create our spatial points object
melesmeles.sp=st_as_sf(melesmeles.latlong,coords=c("x","y"),crs="epsg:4326")

plot(melesmeles.sp)



#First set the extent to the study area
scot=st_read('scotSamp.shp')

#load in the land cover map and then clip to the polygon
LCM=rast("LCMUK.tif")
#crop to the extent of the study area plus a little more (because we will lose a small amount of data in the next step)
LCM=crop(LCM,st_buffer(scot, dist= 1000))
#aggregate LCM raster
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
#project squirrel data
melesmeles.sp=st_transform(melesmeles.sp,crs(LCM))

#now crop our points to the study area
melesmelesFin=melesmeles.sp[scot,]

#finally, mask the LCM to this boundary

LCM=crop(LCM,scot,mask=TRUE)
plot(LCM)
plot(melesmelesFin$geometry,add=T)