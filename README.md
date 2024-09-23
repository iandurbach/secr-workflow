# secr-workflow
A set of scripts for running SCR analyses using the R package `secr`.

Data required:

- List of detections (see **data/example-captfile.csv** for example)
- Detector locations (see **data/example-trapfile.csv** for example)
- Shape file of survey region for extrapolation (see **data/survey-region/Survey_Region.shp** for example)
- Spatial covariates (currently only in raster format, but others possible; see **data/spatial-covs/...** for examples)

Then run code:

* *01_make_capthist.R*: creates input objects for `secr.fit`.
* *02_make_masks.R*: creates an secr mask that covers the whole survey region. Various masks can then be easily constructed from this full mask by subsetting.
* *03_mask_ms_masks.R*: creates multisession secr masks for model fitting.
* *04_model_fitting.R*: fits secr models
* *05_model_results.R*: generate abundances estimates from fitted models and makes some useful plots.

Note: these scripts cover some elements that commonly arise in basic-to-intermediate SCR analyses. They are aimed at countrywide estimates of animal abundance that make use of multiple independent surveys (hence a focus on multisession analyses). They don't aim to cover every possibility. I aim to add some additional functionality to these scripts over time to make them more flexible, but the aim is to keep the workflow simple, rather than comprehensive. Feedback is welcome.