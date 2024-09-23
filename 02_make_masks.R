# Script # 2 / 5
# ###################################################
# The purpose of this script is to create an secr mask that covers the whole survey region.
# Various masks can then be easily constructed from this "full" mask by subsetting.
# This script requires paths to any spatial rasters/shape files to be used as covariates.
# Script to run before this one: 01_make_capthist.R
# Script to run after this one: 03_make_ms_masks.R
# ###################################################

# required libraries, install.packages("package_name") if necessary
library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(terra)
library(secr)
library(sf)
library(sfheaders)
library(tidyr)
library(ggplot2)

select <- dplyr::select
my_crs <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs"
mask_buffer <- 30000

# Extrapolation area: this contains all potential habitat in a country 
extrap_region <- st_read("data/Survey_Region.shp")
extrap_region <- extrap_region %>% st_zm() %>% st_transform("+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs") %>% st_geometry()
plot(extrap_region)

# Survey mask (any points within mask_buffer of a trap, this may contain areas outside the 
# extrapolation area e.g. outside the country if trap close to border and no barriers)

## load traps
load("output/example_secr_inputs_capthist.Rdata")

## buffer around all traps
traps_sf <- st_as_sf(traps, coords = c("x", "y"), crs = my_crs)
full_mask <- traps_sf %>% st_buffer(mask_buffer) %>% st_union() %>% st_geometry()

# Full region is the union of extrapolation region and survey mask

full_region <- st_union(extrap_region, full_mask)

## check that everything looks correct
plot(full_region, border = "white")
plot(full_mask, border = "green", add = TRUE) 
plot(extrap_region, border = "red", add = TRUE) 
plot(traps_sf, add = TRUE, col = "gray")
plot(full_region, border = "blue", add =TRUE, lty = 2, lwd = 2)

# Now turn this into a secr mask by creating grid of points and add covariates at those points

## add another small buffer to avoid any unexpected effects at the boundaries
tmp <- st_make_grid(st_buffer(full_region, 10000), cellsize = c(2000, 2000), what = "centers") 

## crop to full region and convert to points and from there to secr mask
tmp1 <- tmp[full_region]
tmp1 <- tmp1 %>% st_as_sf() %>% sf_to_df() %>% select(x, y)
mask <- read.mask(data = tmp1, spacing = 2000)
rm(tmp, tmp1)

## checks
plot(mask)
attr(mask, "area") * nrow(mask) / 100 # mask area in km2
summary(mask)

## add spatial covariates

### read in
elev <- rast("data/example-elev.tif") # lat-long
rugged <- rast("data/example-tri.tif") # utm
roadden <- rast("data/example-roadden.tif") # utm

### check resolution and projection
elev
rugged
roadden

### project to UTM for use with secr objects
elev <- project(elev, y = my_crs)
rugged <- project(rugged, y = my_crs)
roadden <- project(roadden, y = my_crs)

### change names if desired
names(elev) <- "elev"
names(rugged) <- "rugged"
names(roadden) <- "roadden"

### most of these rasters are at much higher res than our mask spacing.
### generally don't want to know the covariate value JUST at the mask point.
### the more relevant quantity is the covariate value in the vicinity of the mask point i.e. in the mask "cell".
### so aggregate the raster cells so the resolution is approximately the same as the mask resolution.
### to do this work out distance between raster cells and then implied aggregation factor.
### won't be able to get exactly the same res on everything, because aggregation factor is an integer (n cells to combine).

## current res
elev
rugged
roadden

## aggregate any covariates as needed
elev <- aggregate(elev, fact = floor(c(2,2)))
roadden <- aggregate(roadden, fact = floor(c(2,2)))

## current res
elev
rugged
roadden

## check for NAs in covariate
sum(is.na(terra::extract(elev, as.matrix(mask))))
sum(is.na(terra::extract(rugged, as.matrix(mask))))
sum(is.na(terra::extract(roadden, as.matrix(mask))))

## add covariates to mask
mask <- addCovariates(mask, elev)
mask <- addCovariates(mask, rugged)
mask <- addCovariates(mask, roadden)

## replace missing NAs (mean replacement for variables with few NAs)
for(i in names(covariates(mask))){
  covariates(mask)[, i][is.na(covariates(mask)[, i])] <- mean(covariates(mask)[, i], na.rm = TRUE)
}

# scale covariates, first find mean and std devs over full mask
for(i in names(covariates(mask))){
  newi <- paste0("std_",i)
  covariates(mask)[, newi] <- as.numeric(scale(covariates(mask)[, i]))
}
names(covariates(mask))

# crop to altitude range
mask_in <- (covariates(mask)$elev < 5500) & (covariates(mask)$elev > 2500)
mask_altrange <- subset(mask, subset = mask_in)

# crop further to just within country
tmp <- lapply(st_intersects(st_as_sf(mask, coords = c("x", "y"), crs = my_crs), extrap_region), length) %>% unlist()
mask_country <- subset(mask, subset = (tmp == 1))
tmp <- lapply(st_intersects(st_as_sf(mask_altrange, coords = c("x", "y"), crs = my_crs), extrap_region), length) %>% unlist()
mask_country_altrange <- subset(mask_altrange, subset = (tmp == 1))

save(mask, mask_altrange, mask_country, mask_country_altrange, 
     mask_buffer, full_region, extrap_region, 
     file = paste0("output/example_secr_inputs_mask.RData"))

