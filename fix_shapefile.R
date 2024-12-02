#!/opt/R/4.2.3/bin/Rscript

#### LIBS ####
library(raster) # For convenience, shapefile function and show method for Spatial objects
library(rgeos)
library(magrittr)
library(rgdal)

#### VARS ####
layer = "TM_WORLD_BORDERS-0.3" # put the name of the file to correct
shape <- readOGR(dsn=path.expand("shape"),layer=layer) # shape has to be the folder name

#### MAIN ####
gIsValid(shape)
shape <- gBuffer(shape, byid = T, width = 0)

writeOGR(obj = shape, dsn = getwd(), layer = layer, driver = "ESRI Shapefile")
