# Script # 3 / 5
# ###################################################
# The purpose of this script is to fit secr models
# Script to run before this one: 03_make_ms_masks.R
# Script to run after this one: 05_model_results.R
# ###################################################

library(secr)

load("output/example_secr_inputs_mask.RData")
load("output/example_secr_inputs_msmask.RData")
load("output/example_secr_inputs_capthist.Rdata")

names(covariates(mask_ms_altrange)[[1]])
names(covariates(traps_ms)[[1]])
summary(ch_ms, terse = TRUE)

# autoini specifies which session to use to generate starting values. The default (session 1)
# may need to be changed if s1 has no or few detections.

# fit null model
m0 <- secr.fit(ch_ms,
               detectfn="HHN",
               mask = mask_ms_altrange,
               model=list(D ~ 1, lambda0 ~ 1, sigma ~ 1),
               verify = FALSE, details = list(autoini = 1))

# model building for density covariates
m1 <- secr.fit(ch_ms,
               detectfn="HHN",
               mask = mask_ms_altrange,
               model=list(D ~ std_elev, lambda0 ~ 1, sigma ~ 1),
               verify = FALSE,
               start = m0)

m2 <- secr.fit(ch_ms,
                detectfn="HHN",
                mask = mask_ms_altrange,
                model=list(D ~ std_rugged, lambda0 ~ 1, sigma ~ 1),
                verify = FALSE,
                start = m0)

m3 <- secr.fit(ch_ms,
               detectfn="HHN",
               mask = mask_ms_altrange,
               model=list(D ~ std_roadden, lambda0 ~ 1, sigma ~ 1),
               verify = FALSE,
               start = m0)
# check AIC
AIC(m0, m1, m2, m3)

# model building for encounter covariates
# note: can search all combinations of density and encounter covariates, but can
# be time consuming
m1.1 <- secr.fit(ch_ms,
                 detectfn="HHN",
                 mask = mask_ms_altrange,
                 model=list(D ~ std_elev, lambda0 ~ reconyx, sigma ~ 1),
                 verify = FALSE,
                 start = m1)

m1.2 <- secr.fit(ch_ms,
                detectfn="HHN",
                mask = mask_ms_altrange,
                model=list(D ~ std_elev, lambda0 ~ lure_01, sigma ~ 1),
                verify = FALSE,
                start = m1)

m1.3 <- secr.fit(ch_ms,
                 detectfn="HHN",
                 mask = mask_ms_altrange,
                 model=list(D ~ std_elev, lambda0 ~ habitat, sigma ~ 1),
                 verify = FALSE,
                 start = m1)

AIC(m1, m1.1, m1.2, m1.3)

save(list = ls(pattern = "m[0-9]+"), file = paste0("output/example_fitted_models.Rdata"))

