#!/opt/R/4.2.3/bin/Rscript

                            ###Script to make maps###

#"C:\Program Files\R\R-4.1.2\bin\Rscript.exe" Map_making.R --extent (-8,2,50,59) --out Sampling_map_UK.png --coordinates UK_sampling_coordinates.csv

# "C:\Program Files\R\R-4.1.2\bin\Rscript.exe" Map_making.R --extent="(-10, 5.1, 40.5, 59)" --coordinates UK-Cat_sampling_coordinates.csv --out Sampling_map_both.png

# "C:\Program Files\R\R-4.1.2\bin\Rscript.exe" Map_making.R --extent="(-10, 15, 40.5, 59)" --coordinates UK-Cat_sampling_coordinates.csv --out Sampling_map_both.png

#### Dependencies ####
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(ggmap)
library(raster)
library(rgdal)
library(sp)
library(rgeos)
library(Hmisc)
library(optparse)

#### ARGS ####
option_list=list(
  make_option(c("-e", "--extent"), type="character", default=NULL,
              help="Extent of the map to plot, seprated by commas, in this order: longitude lower left corner, longitude top right, latitude lower left, latitude top right",
              metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="Map.png",
              help="Name to use for output files",
              metavar="character"),
  
  make_option(c("-c", "--coordinates"), type="character", default=NULL,
  
              help="Name of the file containing the coordinates (must be semicolon-separated CSV)",
              metavar="character"),
  make_option(c("--color"), type="character", default = "grey",
              help="Name of the column to use to color the points.",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#### OPERATORS ####
`%!in%` = Negate(`%in%`)

#### ARG CHECKS ####
if (!("coordinates" %in% names(opt))){
  print("ERROR: Please provide a coordinates file")
  quit()
}

#### VARS ####
# world_HD <- "shape/gadm41_GBR_0" # Breat Britain borders
world <- "TM_WORLD_BORDERS-0.3"
# world_bad <- "TM_WORLD_BORDERS_SIMPL-0.3"
layer <- world # Select this name from the objects above, based on the quality you prefer. The low-quality one takes less time to be plotted
# adj <- "gadm41_IRL_0" # Ireland borders
# adj2 <- "gadm36_FRA_0" # France borders
# adj3 <- "gadm41_IMN_0" # Isle of Man borders
# adj4 <- "gadm36_ESP_0" # Spain borders
# adj5 <- "Catalonia_shapefile" # Catalonia

if ("extent" %in% names(opt)){
  ext <- gsub(pattern = "(", replacement = "", x = opt$extent, fixed = T)
  ext <- gsub(pattern = ")", replacement = "", x = ext, fixed = T)
  ext <- as.double(as.vector(strsplit(ext, split=",")[[1]]))
}

out <- opt$out
coords <- opt$coordinates
color <- opt$color

#### MAIN ####
### Read dataset ###
dataset <- read.csv(coords, header = TRUE, sep = ";",) #Select your .csv with Longitude column called Lon and Latitude column called Lat, in decimal degrees

if ("extent" %!in% names(opt)){
  lowleftlon <- min(dataset$Lon, na.rm = T)
  toprightlon <- max(dataset$Lon, na.rm = T)
  lowleftlat <- min(dataset$Lat, na.rm = T)
  toprightlat <- max(dataset$Lat, na.rm = T)
  
  ext <- c(lowleftlon, toprightlon, lowleftlat, toprightlat)
}

### plot map ###
shape <- readOGR(dsn=path.expand("shape"),layer=layer)
shape <- crop(shape, extent(ext))
# adjshape<-readOGR(dsn=path.expand("shape"),layer=adj)
# adjshape<-crop(adjshape, extent(ext))
# adj2shape <- readOGR(dsn=path.expand("shape"),layer=adj2)
# adj2shape<-crop(adj2shape, extent(ext))
# adj3shape <- readOGR(dsn=path.expand("shape"),layer=adj3)
# adj3shape<-crop(adj3shape, extent(ext))
# adj4shape <- readOGR(dsn=path.expand("shape"),layer=adj4)
# adj4shape<-crop(adj4shape, extent(ext))
# adj5shape <- readOGR(dsn=path.expand("shape"),layer=adj5)
# adj5shape<-crop(adj5shape, extent(ext))


# pic <- stack("HYP_HR_SR_W.tif")
pic <- getMap(resolution = "high")
pic <- gBuffer(pic, byid = T, width = 0)
pic <- crop(pic, extent(ext))


# pdf(out, height=5, width=5)
png(out, height=2400, width=3100, res=600)
par(mar=c(5,5,2,2)+0.1, cex=0.75, oma=c(3.5, 3, 1, 1))
#mgp.axis.labels(c(1, 0.2, 1), type='y')
#mgp.axis.labels(c(1.0, 0.2, 0), type='x')
# plotRGB(pic, axes=F) +
plot(pic, border = NA, col = "lightgrey") + 
  plot(shape, border="black",lwd=1, add = TRUE)
  # plot(adj5shape, border="black",lwd=1, add = TRUE)
  # plot(adjshape, border="black",lwd=1, add = TRUE) + plot(adj2shape, border="black", lwd=1, add=T) +
  # plot(adj3shape, border = "black", lwd=1, add=T) +
  # plot(adj4shape, border = "black", lwd=1, add=T) +
# mgp.axis(1, at=seq(from=ext[1], to=ext[2], by=(ext[2]-ext[1])/4), labels=round(seq(from=xmin(pic), 
#                                                                                    to=xmax(pic), 
#                                                                                    by = (xmax(pic)-xmin(pic))/4), digits = 1), 
#          axistitle = "Longitude (\u00B0E)", line=0)
# mgp.axis(2, at=seq(from=ext[3], to=ext[4], by=(ext[4]-ext[3])/3), labels=round(seq(from=ymin(pic), 
#                                                                                    to=ymax(pic), 
#                                                                                    by = (ymax(pic)-ymin(pic))/3), digits = 1), 
#          axistitle="Latitude (\u00B0N)", line=0, las=1)
# scalebar(100, xy = c(ext[1]+0.75*(ext[2]-ext[1]), ext[3]+0.9*(ext[4]-ext[3])), type = "bar", below = "kilometers")

### plot records ###
points(dataset$Lon,dataset$Lat,pch=21, cex=4/3, col=color, bg=color) #you can play with colors (coll and bg), shapes (pch) and sizes (cex).

# mtext("Longitude (\u00B0E)", line=2, side=1, cex=1, outer=T)
# mtext("Latitude (\u00B0N)", line=-6.3, side=2, cex=1, las=0, outer=T)

dev.off()

quit()