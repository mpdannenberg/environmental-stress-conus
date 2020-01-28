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






## Notes
*ITRDB data can be downloaded through NOAA’s National Centers for Environmental Information and the World Data Center for Paleoclimatology (https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring)

**DBH measurements are available in the data folder. Ring widths for Washington sites are available from the ITRDB: BCR=wa136, FOD=wa138, RRR=wa139, SND=wa140, SUG=wa141, TMR=wa142.

***Ring widths for site EFV and DBH measurements for all sites are available in the data folder. Ring widths for remaining Wyoming sites are available from the ITRDB: CWC=wy061, ETC=wy063

****Ring widths and DBH measurements for Harvard Forest sites are available from the Harvard Forest Data Archive (https://harvardforest.fas.harvard.edu/harvard-forest-data-archive)

