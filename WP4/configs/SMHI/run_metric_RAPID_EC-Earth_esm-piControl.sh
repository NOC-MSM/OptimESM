#!/usr/bin/env bash

GLOBIGNORE="*"
# ------------------------------------------------------------------------------
# run_metric_RAPID_EC-Earth_esm-piControl.sh
# 
# Description: Script to run METRIC RAPID AMOC AMOC diagnostic for MOHC EC-Earth
# esm-piControl OptimESM simulation.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
#
# Created On: 2025-10-31
#
#
# Note: Use `conda activate env_metric` prior to running this script.
# ------------------------------------------------------------------------------

# -- Input arguments to ./run_metric -- #
# Filepaths to config file:
fpath_config=/g100_work/optim_IAC/research/noc/otooth/OptimESM_WP4/configs/SMHI/config_RAPID_EC-Earth.ini

# Filepaths to eORCA1 monthly mean output files:
exp_id=r1i1p1f1
exp_name=esm-piControl

fdir="/g100_store/DRES_OptimESM/ESGF/prepub/smhi/CMIP6Plus/CMIP/EC-Earth-Consortium/EC-Earth3-ESM-1/${exp_name}/${exp_id}"
fpath_t="${fdir}/OPmonLev/thetao/gn/v*/thetao_OPmonLev_EC-Earth3-ESM-1_${exp_name}_${exp_id}_gn_*.nc"
fpath_s="${fdir}/OPmonLev/so/gn/v*/so_OPmonLev_EC-Earth3-ESM-1_${exp_name}_${exp_id}_gn_*.nc"
fpath_v="${fdir}/OPmonLev/vo/gn/v*/vo_OPmonLev_EC-Earth3-ESM-1_${exp_name}_${exp_id}_gn_*.nc"
fpath_ssh="${fdir}/OPmon/zos/gn/v*/zos_OPmon_EC-Earth3-ESM-1_${exp_name}_${exp_id}_gn_*.nc"
fpath_taux="${fdir}/OPmon/tauuo/gn/v*/tauuo_OPmon_EC-Earth3-ESM-1_${exp_name}_${exp_id}_gn_*.nc"

# Update the output filename in the config .ini file:
outfname="EC-Earth3-ESM-1_${exp_name}_${exp_id}"
sed -i "s/^name *= *.*/name = ${outfname}/" $fpath_config

# -- Run METRIC -- #
echo "In Progress: Calculating RAPID 26N AMOC Diagnostics..."

metric run -c $fpath_config -t $fpath_t -s $fpath_s -v $fpath_v -ssh $fpath_ssh -taux $fpath_taux

echo "Completed: Calculated RAPID 26N AMOC Diagnostics."
