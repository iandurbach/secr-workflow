# Script # 1 / 5
# ###################################################
# The purpose of this script is to create the input objects for secr.fit.
# This consists of a "capthist" object containing capture histories and trap locations.
# This script takes one csv file for detector locations and one csv file for detections and creates the secr "capthist" object.
# Script to run before this one: none.
# Script to run after this one: 02_make_masks.R
# ###################################################

# required libraries, install.packages("package_name") if necessary
library(readxl)
library(dplyr)
library(secr)
library(tidyr)
library(sf)
library(sfheaders)
library(stringr)

# traps 

## read in
traps <- read.csv("data/example-trapfile.csv")

# encounters

## read in
capts <- read.csv("data/example-captfile.csv")

# make (multisession) secr inputs

## traps object (all sessions, including empty ones)

traps_list <- split(data.frame(traps), f = traps$Session)  
names(traps_list)
traps_ms <- list()
for(i in 1:length(traps_list)){
  traps_covs <- traps_list[[i]] %>% dplyr::select(camera_type, lure_01, habitat, terrain)
  traps_basic <- traps_list[[i]] %>% dplyr::select(ID, x, y, Effort)
  traps_ms[[i]] <- read.traps(data = traps_basic, detector = "count", trapID = "ID", binary.usage = FALSE)
  covariates(traps_ms[[i]])$camera_type <- factor(traps_covs$camera_type, levels = unique(traps$camera_type))
  covariates(traps_ms[[i]])$reconyx <- ifelse(covariates(traps_ms[[i]])$camera_type == "reconyx", 1, 0)
  covariates(traps_ms[[i]])$lure_01 <- factor(traps_covs$lure_01, levels = unique(traps$lure_01))
  covariates(traps_ms[[i]])$habitat <- factor(traps_covs$habitat, levels = unique(traps$habitat))
  covariates(traps_ms[[i]])$terrain <- factor(traps_covs$terrain, levels = unique(traps$terrain))
}
head(covariates(traps_ms[[1]]))
attr(traps_ms[[1]], "usage")
traps_ms <- shareFactorLevels(traps_ms, stringsAsFactors = FALSE)

## capthist object 

### add empty sessions
emptysess <- setdiff(unique(traps$Session), unique(capts$Session))
emp_df <- data.frame(Session = emptysess, ID = "NONE", occasion = 1, trapID = 0)
capts <- rbind(capts, emp_df)
ch_ms <- make.capthist(captures = capts, traps = traps_ms)
verify(ch_ms)
summary(ch_ms, terse = TRUE)

save(capts, traps, traps_ms, ch_ms,
     file = paste0("output/example_secr_inputs_capthist.Rdata"))


