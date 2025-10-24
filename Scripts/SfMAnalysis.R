#### SfM Photogrammetry Analysis #####
### Created by: Ashley Carreiro
### Updated on: 2025_09-24
#update

#Load Libraries ####
library(tidyverse)
library(here)
library(devtools)
library(habtools)
library(raster)
library(sf)

#### Load in DEMS ####

Oki11_Dem <- raster("Data/SfM/Oki11_DEM.tif")
plot(Oki11_Dem)

#### Take a subset DEM of size = 2 ####

oki11_dem_Square <- dem_crop(Oki11_Dem, x0 = -6, y0 = 7, L = 0.5, plot =TRUE)

#Take the surface area of the square
surface_area(oki11_dem_Square)


#Get the rugosity of the square
rg(oki11_dem_Square)

#Height Range
hr(oki11_dem_Square)

#Fractal Dimension
fd(oki11_dem_Square, method = "sd", lvec = c(0.31, 0.063, 0.125, 0.25, 0.5), diagnose = TRUE, parallel = TRUE)

#DEM Split 

dem_list <- dem_split(oki11_dem_Square, 0.1)
length(dem_list)
plot(Oki11_Dem)
rect(0, 0, 0.5, 0.5)
plot(dem_list[[1]], add=TRUE, legend = FALSE) #don't know hwat this is doing 
plot(dem_list[[25]], add = TRUE, legend = FALSE)


#Metrics for dem list

rdhs <- lapply(dem_list, rdh, method_fd = "sd", lvec = c(0.031, 0.063, 0.125, 0.25, 0.5), parallel = TRUE)
rdhs <- rdhs %>%
  plyr::ldply()


#Plot RDH 

ggplot(data = rdhs, aes(x = R, y = H, color = D, size = D)) +
  geom_point()


#Sample DEM 
dem_sample(oki11_dem_Square, L = 0.2, plot = TRUE)


