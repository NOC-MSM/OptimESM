#!/usr/bin/env bash

GLOBIGNORE="*"
# ------------------------------------------------------------------------------
# run_metric_RAPID_EC-Earth_historical.sh
# 
# Description: Script to run METRIC RAPID AMOC AMOC diagnostic for MOHC EC-Earth
# historical OptimESM simulation.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
#
# Created On: 2025-03-01
#
#
# Note: Use `conda activate env_metric` prior to running this script.
# ------------------------------------------------------------------------------

# -- Input arguments to ./run_metric -- #
# Filepaths to config file:
fpath_config=/g100_work/optim_IAC/research/noc/otooth/OptimESM_WP4/configs/IPSL/config_RAPID_IPSL_esm-hist.ini

# Filepaths to eORCA1 monthly mean output files:
exp_id=r1i1p3f1
fdir="/g100_store/DRES_OptimESM/ESGF/external/20250619/CMIP6Plus/CMIP/IPSL/IPSL-CM6-ESMCO2/esm-hist/${exp_id}/Omon"
fpath_t="${fdir}/thetao/gn/v*/thetao_Omon_IPSL-CM6-ESMCO2_esm-hist_${exp_id}_gn_*.nc"
fpath_s="${fdir}/so/gn/v*/so_Omon_IPSL-CM6-ESMCO2_esm-hist_${exp_id}_gn_*.nc"
fpath_v="${fdir}/vo/gn/v*/vo_Omon_IPSL-CM6-ESMCO2_esm-hist_${exp_id}_gn_*.nc"
fpath_ssh="${fdir}/zos/gn/v*/zos_Omon_IPSL-CM6-ESMCO2_esm-hist_${exp_id}_gn_*.nc"
fpath_taux="${fdir}/tauuo/gn/v*/tauuo_Omon_IPSL-CM6-ESMCO2_esm-hist_${exp_id}_gn_*.nc"

# -- Run METRIC -- #
echo "In Progress: Calculating RAPID 26N AMOC Diagnostics..."

metric run -c $fpath_config -t $fpath_t -s $fpath_s -v $fpath_v -ssh $fpath_ssh -taux $fpath_taux

echo "Completed: Calculated RAPID 26N AMOC Diagnostics."
