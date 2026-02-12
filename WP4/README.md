# OptimESM WP4 - RAPID 26.5N Atlantic Meridional Overturning Circulation Diagnostics

## Description:
Observational-equivalent Atlantic Meridional Overturning Circulation (AMOC) diagnostics calculated using the latest version [2025-08-18] of the Meridional ovErTurning ciRculation diagnostIC (METRIC) package (see Danabasoglu, G., Castruccio, F. S., Small, R. J., Tomas, R., Frajka-Williams, E., and Lankhorst, M., (2021). Revisiting AMOC Transport Estimates from Observations and Models. Geophysical Research Letters).

The METRIC package calculates the AMOC & its contributing flow components (Florida Current, WBW, interior, Ekman etc.), Meridional Heat & Freshwater and property sections as done in RAPID-MOCHA observations. Model "truth" diagnostics calculated by integrating the meridional velocity field over longitude and depth dimensions are also included, although it should be noted this is done by approximation rather than using vo(t, lev, i) * e1v(i) * e3v(t, lev, i) for j=j_rapid_26N as formulated in NEMO.

This directory contains:
* `metric/` - metric Python package cloned from GitHub (git clone git@github.com:oj-tooth/metric.git)
* `configs/` - directory of metric RAPID 26.5N configuration `.ini` files, run scripts and pre-processing Jupyter Notebooks (this include definitions of RAPID 26.5N meridional sections in the native model i,j coordinates and t/u/v mask file creation) divided into separate institutions.
* `data/` - metric output netCDF files for the historical simulations (1850-2014) available on CINECA.
* `analysis/` - Jupyter Notebooks perfoming a validation of OptimESM 2004-2022 RAPID 26.5N diagnostics with observations.

**RAPID 26.5N AMOC diagnostics have been performed for the following historical simulations:**

### **UKESM1**
*/g100_store/DRES_OptimESM/ESGF/prepub/mohc/20240619/CMIP6/CMIP/MOHC/UKESM1-2/esm-hist/r1i1p1f1/Omon/
*/g100_store/DRES_OptimESM/ESGF/prepub/mohc/20241218/CMIP6/CMIP/MOHC/UKESM1-2/esm-hist/r2i1p1f1/Omon/
*/g100_store/DRES_OptimESM/ESGF/prepub/mohc/20241218/CMIP6/CMIP/MOHC/UKESM1-2/esm-hist/r3i1p1f1/Omon/

### **EC-Earth**
* /g100_store/DRES_OptimESM/ESGF/prepub/cnr/20241129/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-ESM-1/esm-hist/r1i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/prepub/bsc/20250212/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-ESM-1/esm-hist/r2i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/prepub/dmi/20240814/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-ESM-1/esm-hist/r3i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/prepub/ulund/hist_r4/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-ESM-1/esm-hist/r4i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/prepub/smhi/20240622/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-ESM-1/esm-hist/r5i1p1f1/Omon/

### **IPSL**
* /g100_store/DRES_OptimESM/ESGF/external/20250619/CMIP6Plus/CMIP/IPSL/IPSL-CM6-ESMCO2/esm-hist/r1i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/external/20250619/CMIP6Plus/CMIP/IPSL/IPSL-CM6-ESMCO2/esm-hist/r2i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/external/20250619/CMIP6Plus/CMIP/IPSL/IPSL-CM6-ESMCO2/esm-hist/r3i1p1f1/Omon/
* /g100_store/DRES_OptimESM/ESGF/external/20250619/CMIP6Plus/CMIP/IPSL/IPSL-CM6-ESMCO2/esm-hist/r4i1p1f1/Omon/

### **CNRM**
No metric RAPID 26.5N AMOC diagnostics currently available due to missing variables:
* /g100_store/DRES_OptimESM/ESGF/prepub/cnrm/20240805/CMIP6/CMIP/CNRM-CERFACS/CNRM-ESM2-1/esm-hist/r1i1p2f2/Omon