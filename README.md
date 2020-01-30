# environmental-stress-conus
Delineating environmental stresses to primary production of U.S. forests from tree rings: Effects of climate seasonality, soil, and topography


## Calculate stress index and conduct sensitivity tests

### Base functions for stress index
StressIndexFunctions.R


### Calculate stress index for ITRDB sites*
calculate_stress_index_itrdb.R


### Conduct sensitivity analyses using test sites in Washington**, Wyoming***, and Harvard Forest**** 
test_sensitivity_harvard.R

test_sensitivity_pipo_psme.R

make_sensitivity_figs.m


## Read and process stress index and environmental data into MATLAB structures
### Calculate potential evapotranspiration (PET) from PRISM data*****
calc_prism_pet.m


### Read tree-ring stress index and environmental data into MATLAB structure.
treesi_read.m

treesi_read_environmental_data.m


## Read and process relative stress index from ITRDB data
### Read ITRDB metadata and detrend with spline
read_itrdb.r

### Read detrended ITRDB chronologies into MATLAB structure
read_itrdb.m

### Read environmental data into ITRDB MATLAB structure
treesi_read_environmental_data_sr.m


## Climate, soil, and topography analyses
### Calculate correlation coefficients and calibrate RF models for each ecoregion and each variable subset (WARNING: very slow and create new folders and lots of files)
calibrate_rf_sstar_C.m

calibrate_rf_sstar_SC.m

calibrate_rf_sstar_TC.m

calibrate_rf_sstar_TSC.m


### Calculate correlation coefficients and calibrate RF models for each ecoregion and each variable subset FOR THE RELATIVE STRESS INDEX (Sr) (WARNING: very slow and create new folders and lots of files)
calibrate_rf_sr_C.m

calibrate_rf_sr_SC.m

calibrate_rf_sr_TC.m

calibrate_rf_sr_TSC.m


### Fit linear mixed effects models with interactions
calibrate_lme.m


### Make correlation figures
make_correlation_maps.m


### Make figures and tables of random forest diagnostics and variable importance
make_rf_model_skill_map.m

make_rf_variable_importance_figure.m

build_table_rf_sstar.m


## Miscellaneous figures
### Map of ITRDB sites used in the analysis
make_itrdb_sites_map.m


### Map of sample depth by site
make_supplemental_sample_depth_map.m


## Notes
*ITRDB data can be downloaded through NOAA’s National Centers for Environmental Information and the World Data Center for Paleoclimatology (https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring)

**DBH measurements are available in the data folder. Ring widths for Washington sites are available from the ITRDB: BCR=wa136, FOD=wa138, RRR=wa139, SND=wa140, SUG=wa141, TMR=wa142.

***Ring widths for site EFV and DBH measurements for all sites are available in the data folder. Ring widths for remaining Wyoming sites are available from the ITRDB: CWC=wy061, ETC=wy063

****Ring widths and DBH measurements for Harvard Forest sites are available from the Harvard Forest Data Archive (https://harvardforest.fas.harvard.edu/harvard-forest-data-archive)

*****Monthly PRISM climate data are available from http://www.prism.oregonstate.edu/. Code for FAO Penman-Monteith method is available at: https://github.com/mpdannenberg/matts-matlab-code/blob/master/fao_pm.m.

