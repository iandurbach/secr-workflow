# Script # 5 / 5
# ###################################################
# The purpose of this script is to generate abundance estimates from fitted models and make some density plots
# Script to run before this one: 04_model_fitting.R
# Script to run after this one: none
# ###################################################

library(dplyr)
library(secr)
library(sf)
library(ggspatial)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)
library(rlang)

select <- dplyr::select
my_crs <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs"

# Basemap for plots
esri_topo <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                    'World_Topo_Map/MapServer/tile/${z}/${y}/${x}.jpeg')

# load output from previous scripts
load("output/example_secr_inputs_mask.Rdata")
load("output/example_secr_inputs_msmask.RData")
load("output/example_secr_inputs_capthist.Rdata")
load("output/example_fitted_models.Rdata")

nsess <- ncol(summary(m0$capthist, terse = TRUE))
sessions_included <- colnames(summary(m0$capthist, terse = TRUE))
ch_ms <- subset(ch_ms, sessions = sessions_included)
mask_ms_altrange <- mask_ms_altrange[sessions_included]

cellarea_km2 <- attr(mask_country_altrange, "area") / 100
cellarea_km2

#########################################
### Table of abundance, density, AICs ###
#########################################

# abundance and density over all country habitat

## specify which mask defines "all country habitat"
extrap_mask <- mask_country_altrange

## confirm elevation limits
summary(covariates(extrap_mask)[,"elev"])

## abundance
allmods <- ls(pattern = "^m[0-9]+")
R.N.slmod <- list()
for(i in allmods){
  mod <- sym(i)
  R.N.slmod[[i]] <- data.frame(modelname = i, region.N(eval(mod), region = extrap_mask)[[1]][1,1:4])
}
R.N.slmods <- do.call(rbind, R.N.slmod)

## mean density = abundance / area
R.N.slmods$D <- R.N.slmods$estimate / (nrow(extrap_mask) * cellarea_km2 / 100)
R.N.slmods$Dlcl <- R.N.slmods$lcl / (nrow(extrap_mask) * cellarea_km2 / 100)
R.N.slmods$Ducl <- R.N.slmods$ucl / (nrow(extrap_mask) * cellarea_km2 / 100)
R.N.slmods$CV <- round(R.N.slmods$SE.estimate / R.N.slmods$estimate * 100, 0)

# abundance and density over mask region only

## add variable specifying if extrapolation mask point is within a mask used to fit model
extrap_mask$in_sess_masks <- lapply(st_intersects(st_as_sf(extrap_mask, coords = c("x", "y"), crs = my_crs), all_sess_masks), length) %>% unlist()
in_sess_mask <- subset(extrap_mask, subset = extrap_mask$in_sess_masks == 1)

## abundance
allmods <- ls(pattern = "^m[0-9]+")
R.N.slmod <- list()
for(i in allmods){
  mod <- sym(i)
  R.N.slmod[[i]] <- data.frame(modelname = i, region.N(eval(mod), region = in_sess_mask)[[1]][1,1:4])
}
R.N.slmods_mask <- do.call(rbind, R.N.slmod)

## mean density = abundance / area
R.N.slmods_mask$D <- R.N.slmods_mask$estimate / (nrow(in_sess_mask) * 4/100)
R.N.slmods_mask$Dlcl <- R.N.slmods_mask$lcl / (nrow(in_sess_mask) * 4/100)
R.N.slmods_mask$Ducl <- R.N.slmods_mask$ucl / (nrow(in_sess_mask) * 4/100)
R.N.slmods_mask$CV <- round(R.N.slmods_mask$SE.estimate / R.N.slmods_mask$estimate * 100, 0)

## rename to distinguish between full and reduced masks
names(R.N.slmods_mask) <- paste0(names(R.N.slmods_mask),"_msk")
R.N.slmods <- cbind(R.N.slmods, R.N.slmods_mask)

# add AICs
AIC.slmods <-  AIC(m0, m1, m1.1, m1.2, m1.3, m2, m3)
AIC.slmods$modelname = row.names(AIC.slmods)

# combine results and save
AIC.slmods <- AIC.slmods %>% left_join(R.N.slmods, by = "modelname") 
write.csv(AIC.slmods %>% dplyr::select(modelname, model, AICc, dAICc, AICcwt, estimate, lcl, ucl, D, Dlcl, Ducl, CV,
                                       estimate_msk, lcl_msk, ucl_msk, D_msk, Dlcl_msk, Ducl_msk, CV_msk), 
          file = paste0("output/results/example_results_abundance_AIC.csv"))

#############################################################
### Map of the survey region and expected SL AC densities ###
#############################################################

# turn traps into an sf object
traps_sf <- st_as_sf(traps, coords = c("x", "y"), crs = my_crs)

# map of survey region
g0 <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 8) +
  layer_spatial(st_geometry(all_sess_masks), colour = "red", lwd = 0.5, fill = NA) + 
  layer_spatial(st_geometry(extrap_region), colour = "purple", lwd = 0.5, fill = NA) + 
  layer_spatial(traps_sf, size = 0.25, col = "black", alpha = 0.5) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")

g0

# Extract density estimates from user-specified model
bestmod <- m1
tmp <- predictDsurface(bestmod, mask = extrap_mask, parameter = "D", cl.D = FALSE)

# Make data frame with density, using a grid of cells for nicer plotting
cellSize <- attr(extrap_mask, "spacing")
exp_ac_dens <- st_as_sf(data.frame(extrap_mask), coords = c("x", "y"), crs = my_crs) 
exp_ac_dens <- st_buffer(exp_ac_dens, dist = cellSize/2, endCapStyle="SQUARE", nQuadSegs = 1)
exp_ac_dens$D_per_100km2_v1 <- covariates(tmp)[,"D.0"] * 10000

# Plot expected activity center density
g1_c <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 8) +
  layer_spatial(exp_ac_dens %>% 
                  dplyr::filter(!is.na(D_per_100km2_v1)), 
                aes(colour = D_per_100km2_v1, fill = D_per_100km2_v1), size = 0, alpha = 0.8) + 
  layer_spatial(st_geometry(extrap_region), colour = "red", lwd = 0.5, fill = NA) + 
  layer_spatial(st_geometry(all_sess_masks), colour = "white", lwd = 0.35, lty = 2, fill = NA) + 
  layer_spatial(traps_sf, col = "red", size = 0.01, alpha = 0.2) +
  scale_fill_viridis_c(option = "D", name = bquote("Density / 100km"^2)) +
  scale_colour_viridis_c(option = "D", guide = "none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")

g1_c

# Plot discretised expected AC density categories

## choose discretisation here (SL/100km2)
quantile(exp_ac_dens$D_per_100km2_v1, probs = seq(0,1,0.1))
my_breaks <- c(0, 0.1, 0.3, 0.75, 1.5, 100)

## making labels
n_breaks <- length(my_breaks)
my_breaks_labels <- c(paste0("<",my_breaks[2]), paste0(my_breaks[2:(n_breaks-2)], rep("-", n_breaks - 3), my_breaks[3:(n_breaks-1)]), paste0(">",my_breaks[n_breaks-1]))

## discretise density
exp_ac_dens_cat <- exp_ac_dens %>% 
  mutate_at(vars(starts_with("D_per_100")), cut, breaks = my_breaks, labels = my_breaks_labels)

## check frequencies in each level
table(exp_ac_dens_cat$D_per_100km2_v1)

## plot
g1_d <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 8) +
  layer_spatial(exp_ac_dens_cat %>% 
                  dplyr::filter(!is.na(D_per_100km2_v1)), 
                aes(colour = D_per_100km2_v1, fill = D_per_100km2_v1), size = 0, alpha = 0.8) + 
  layer_spatial(st_geometry(extrap_region), colour = "red", lwd = 0.5, fill = NA) + 
  layer_spatial(st_geometry(all_sess_masks), colour = "white", lwd = 0.35, lty = 2, fill = NA) + 
  layer_spatial(traps_sf, col = "red", size = 0.01, alpha = 0.2) +
  scale_fill_viridis_d(option = "D", name = bquote("Density / 100km"^2)) +
  scale_colour_viridis_d(option = "D", guide = "none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")

g1_d

## save maps
ggsave("output/results/example_survey_region.png", g0, width = 8, height = 5, dpi = 150)
ggsave("output/results/example_predicted_densities.png", g1_c, width = 8, height = 5, dpi = 150)
ggsave("output/results/example_predicted_densities_cat.png", g1_d, width = 8, height = 5, dpi = 150)

## write shape file with results
tmp <- exp_ac_dens %>% dplyr::select(D = D_per_100km2_v1)
st_write(tmp, "output/results/example_predicted_densities.SHP", delete_dsn = TRUE)

######################################
### Plot density covariate effects ###
######################################

# data prep
names(covariates(extrap_mask))

# plotting function
plot_denscov = function(predmesh, fitmesh, covname, covlabel, logx = FALSE, ymax = 8, xscale = 1){
  
  x <- covariates(predmesh)[, covname]
  meanx <- mean(covariates(predmesh)[, str_remove_all(covname, "std_")])
  sdx <- sd(covariates(predmesh)[, str_remove_all(covname, "std_")])
  preds <- data.frame(id = 1:500)
  preds[, covname] <- seq(min(x), max(x), length.out = 500)
  preds$xt <- seq(min(x), max(x), length.out = 500) * sdx + meanx
  
  predname <- names(preds)[2]
  tmp <- switch(predname,
                std_elev = {tmp <- predict(m1, newdata = preds)},
                std_rugged = {tmp <- predict(m2, newdata = preds)},
                std_roadden = {tmp <- predict(m3, newdata = preds)},
                stop("Choose a different covariate")
  )
  
  preds$est <- 10000 * sapply(tmp, "[", "D","estimate")
  preds$lcl <- 10000 * sapply(tmp, "[", "D","lcl")
  preds$ucl <- 10000 * sapply(tmp, "[", "D","ucl")
  
  x_surveyed <- data.frame(xt = c(covariates(fitmesh)[, covname] * sdx + meanx))
  x_full <- data.frame(xt = c(covariates(predmesh)[, covname] * sdx + meanx))
  
  if(logx){
    p1 <- preds %>% 
      ggplot(aes(x = xt / xscale, y = est)) + 
      geom_line() +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "gray80", alpha = 0.5) + 
      geom_rug(data = x_full, inherit.aes = FALSE, aes(x = xt / xscale), sides = "b", alpha = 0.5, size = 0.25) +
      geom_rug(data = x_surveyed, inherit.aes = FALSE, aes(x = xt / xscale), sides = "t", alpha = 0.5, size = 0.25) +
      coord_cartesian(ylim = c(0, ymax)) +
      scale_x_log10() +
      labs(x = covlabel, y = bquote("Density / 100km"^2)) +
      theme_bw(base_size = 14) + 
      theme(panel.grid = element_blank()) 
  } else {
    p1 <- preds %>% 
      ggplot(aes(x = xt / xscale, y = est)) + 
      geom_line() +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "gray80", alpha = 0.5) + 
      geom_rug(data = x_full, inherit.aes = FALSE, aes(x = xt / xscale), sides = "b", alpha = 0.5, size = 0.25) +
      geom_rug(data = x_surveyed, inherit.aes = FALSE, aes(x = xt / xscale), sides = "t", alpha = 0.5, size = 0.25) +
      coord_cartesian(ylim = c(0, ymax)) +
      labs(x = covlabel, y = bquote("Density / 100km"^2)) +
      theme_bw(base_size = 14) + 
      theme(panel.grid = element_blank()) 
  }
  
  return(p1)
  
}

denscov_elev <- plot_denscov(predmesh = extrap_mask, fitmesh = mask_altrange, 
                             covname = "std_elev", covlabel = "Elevation (km)", ymax = 4, xscale = 1000)

denscov_rugged <- plot_denscov(predmesh = extrap_mask, fitmesh = mask_altrange, 
                             covname = "std_rugged", covlabel = "Ruggedness", ymax = 4, xscale = 1)

denscov_roadden <- plot_denscov(predmesh = extrap_mask, fitmesh = mask_altrange, 
                               covname = "std_roadden", covlabel = "Road density", ymax = 4, xscale = 1)

p2 <- denscov_elev + denscov_rugged + denscov_roadden
p2

## save 
ggsave("output/results/example_density_covs.png", p2, width = 8, height = 3, dpi = 150)
