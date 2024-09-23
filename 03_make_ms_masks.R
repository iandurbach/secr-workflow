# Script # 3 / 5
# ###################################################
# The purpose of this script is to create multisession secr masks for model fitting.
# This is done by cropping full region mask to polygon around traps in each session
# Covariates have already been added, this is just a subsetting
# Script to run before this one: 02_make_masks.R
# Script to run after this one: 04_model_fitting.R
# ###################################################

library(dplyr)
library(raster)
library(secr)
library(tidyr)
library(lubridate)
library(sf)
library(sfheaders)
library(stringr)

select <- dplyr::select
my_crs <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs"

load("output/example_secr_inputs_mask.RData")
load("output/example_secr_inputs_capthist.RData") 
ls() ## what these files contain

## note: session-specific masks are created by cropping the full mask. Each session's mask should be within
## mask_buffer of the traps for that session.

##identify everywhere close to traps in each session
this_mask_poly <- list()
for(i in 1:length(traps_ms)){
  # mask polygon for each session
  this_mask_poly[[i]] <- st_as_sf(traps_ms[[i]], coords = c("x", "y"), crs = my_crs) %>% st_buffer(mask_buffer) %>% st_union()
}

## crop mask points in each session to the polygon just created
mask_sf <- st_as_sf(mask, coords = c("x", "y"), crs = my_crs) %>% st_geometry()
mask_ms <- list()
for(i in 1:length(this_mask_poly)){
  # binary indicator of whether mask points are in mask polygon
  mask_in <- lapply(st_intersects(mask_sf, this_mask_poly[[i]]), length) %>% unlist() %>% as.logical()
  mask_ms[[i]] <- subset(mask, mask_in)
}
class(mask_ms)
class(mask_ms[[1]])
names(covariates(mask_ms[[1]]))
names(mask_ms) <- names(ch_ms)

# crop mask by altitude
mask_ms_altrange <- mask_ms
for(i in 1:length(mask_ms_altrange)){
  mask_in <- covariates(mask_ms_altrange[[i]])$elev < 5500
  mask_ms_altrange[[i]] <- subset(mask_ms_altrange[[i]], subset = mask_in)
}

# create a mask that is the union of all session masks

## union of all session masks
all_sess_masks <- this_mask_poly[[1]]
for(i in 2:length(this_mask_poly)){
  all_sess_masks <- st_union(all_sess_masks, this_mask_poly[[i]])
}

## crop full mask to union of all session masks
tmp <- lapply(st_intersects(st_as_sf(mask, coords = c("x", "y"), crs = my_crs), all_sess_masks), length) %>% unlist()
mask_all_sess <- subset(mask, subset = (tmp == 1))
tmp <- lapply(st_intersects(st_as_sf(mask_altrange, coords = c("x", "y"), crs = my_crs), all_sess_masks), length) %>% unlist()
mask_all_sess_altrange <- subset(mask_altrange, subset = (tmp == 1))

# crop masks just created further to just country
tmp <- lapply(st_intersects(st_as_sf(mask_all_sess, coords = c("x", "y"), crs = my_crs), extrap_region), length) %>% unlist()
mask_all_sess_country <- subset(mask_all_sess, subset = (tmp == 1))
tmp <- lapply(st_intersects(st_as_sf(mask_all_sess_altrange, coords = c("x", "y"), crs = my_crs), extrap_region), length) %>% unlist()
mask_all_sess_country_altrange <- subset(mask_all_sess_altrange, subset = (tmp == 1))

# save list of session mask (sf object)
mask_per_session_sf <- this_mask_poly

save(mask_ms, mask_ms_altrange, 
     mask_all_sess_country, mask_all_sess_country_altrange, 
     mask_per_session_sf, all_sess_masks, 
     file = "output/example_secr_inputs_msmask.Rdata")
