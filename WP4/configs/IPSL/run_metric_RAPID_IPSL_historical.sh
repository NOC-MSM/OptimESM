#!/usr/bin/env bash

GLOBIGNORE="*"
# ------------------------------------------------------------------------------
# run_metric_RAPID_IPSL_historical.sh
# 
# Description: Script to run METRIC RAPID AMOC AMOC diagnostic for IPSL
# IPSL-CM6-ESMCO2 OptimESM simulation.
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
fpath_config=/g100_work/optim_IAC/research/noc/otooth/OptimESM_WP4/configs/IPSL/config_RAPID_IPSL.ini

# Define model & experiment properties:
model_name=IPSL-CM6-ESMCO2
exp_id=r1i1p3f1
exp_name=esm-hist

fdir="/g100_store/DRES_OptimESM/ESGF/external/20250619/CMIP6Plus/CMIP/IPSL/${model_name}/${exp_name}/${exp_id}/Omon"
fpath_t="${fdir}/thetao/gn/v*/thetao_Omon_${model_name}_${exp_name}_${exp_id}_gn_*.nc"
fpath_s="${fdir}/so/gn/v*/so_Omon_${model_name}_${exp_name}_${exp_id}_gn_*.nc"
fpath_v="${fdir}/vo/gn/v*/vo_Omon_${model_name}_${exp_name}_${exp_id}_gn_*.nc"
fpath_ssh="${fdir}/zos/gn/v*/zos_Omon_${model_name}_${exp_name}_${exp_id}_gn_*.nc"
fpath_taux="${fdir}/tauuo/gn/v*/tauuo_Omon_${model_name}_${exp_name}_${exp_id}_gn_*.nc"

# Update the output filename in the config .ini file:
outfname="${model_name}_${exp_name}_${exp_id}"
sed -i "s/^name *= *.*/name = ${outfname}/" $fpath_config

# -- Run METRIC -- #
echo "In Progress: Calculating RAPID 26N AMOC Diagnostics..."

metric run -c $fpath_config -t $fpath_t -s $fpath_s -v $fpath_v -ssh $fpath_ssh -taux $fpath_taux

echo "Completed: Calculated RAPID 26N AMOC Diagnostics."
