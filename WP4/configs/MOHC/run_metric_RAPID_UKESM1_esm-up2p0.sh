#!/usr/bin/env bash

GLOBIGNORE="*"
# ------------------------------------------------------------------------------
# run_metric_RAPID_UKESM1_esm-up2p0-gwl4p0.sh
# 
# Description: Script to run METRIC RAPID AMOC AMOC diagnostic for MOHC UKESM1
# esm-up2p0-gwl4p0 OptimESM simulation.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
#
# Created On: 2025-10-30
#
#
# NOTE: Use `conda activate env_metric` prior to running this script.
# ------------------------------------------------------------------------------

# -- Input arguments to ./run_metric -- #
# Filepaths to config file:
fpath_config="/g100_work/optim_IAC/research/noc/otooth/OptimESM_WP4/configs/MOHC/config_RAPID_UKESM1.ini"

# Filepaths to eORCA1 monthly mean output files:
exp_id=r1i1p1f1
exp_name=esm-up2p0

# -- Defining filepaths & updating config .ini -- #
fdir="/g100_store/DRES_OptimESM/ESGF/prepub/mohc/20250327/CMIP6/CMIP/MOHC/UKESM1-2/${exp_name}/${exp_id}/Omon"
fpath_t="${fdir}/thetao/gn/v*/thetao_Omon_UKESM1-2-LL_${exp_name}_${exp_id}_gn_*.nc"
fpath_s="${fdir}/so/gn/v*/so_Omon_UKESM1-2-LL_${exp_name}_${exp_id}_gn_*.nc"
fpath_v="${fdir}/vo/gn/v*/vo_Omon_UKESM1-2-LL_${exp_name}_${exp_id}_gn_*.nc"
fpath_ssh="${fdir}/zos/gn/v*/zos_Omon_UKESM1-2-LL_${exp_name}_${exp_id}_gn_*.nc"
fpath_taux="${fdir}/tauuo/gn/v*/tauuo_Omon_UKESM1-2-LL_${exp_name}_${exp_id}_gn_*.nc"

# Update the output filename in the config .ini file:
outfname="UKESM1-2-LL_${exp_name}_${exp_id}"
sed -i "s/^name *= *.*/name = ${outfname}/" $fpath_config

# -- Run METRIC -- #
echo "In Progress: Calculating RAPID 26N AMOC Diagnostics..."

metric run -c $fpath_config -t $fpath_t -s $fpath_s -v $fpath_v -ssh $fpath_ssh -taux $fpath_taux

echo "Completed: Calculated RAPID 26N AMOC Diagnostics."
