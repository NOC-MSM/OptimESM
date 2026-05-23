#!/bin/bash

set -euo pipefail
# ================================================================
# run_UKESSM1-2-LL_pipeline.sh
#
# Description: Run UKESM1-2-LL NEMO Pipeline in current process.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-11-11   
# ================================================================

# -- Input arguments to NEMO Pipeline -- #
# Define filepaths:
config_file=osnap/config_UKESM1_esm-hist.toml
log_file=UKESM1-2-LL_osnap_pipeline.log

# Run multiple pipelines:
l_multi=true

# Define Experiment IDs [l_multi=true] -> esm-hist ensemble members...
member_ids=("r1i1p1f1" "r2i1p1f1" "r3i1p1f1")

# -- Python Environment -- #
# Run this script in the env_optimesm conda virtual environment.

if [ "$l_multi" = false ]; then
    # -- Run NEMO Pipeline CLI -- #
    nemo_pipeline run "$config_file" --log "$log_file"

else
    # Iterate over all realisation IDs:
    for member_id in "${member_ids[@]}"; do
        echo "Running ==> $member_id"

        # -- Updating member IDs in config.toml -- #
        sed -i "s|r[0-9]\+i1p1f1|$member_id|g" "$config_file"

        # -- Run NEMO Pipeline CLI -- #
        nemo_pipeline run "$config_file" --log "$log_file"

        echo "Completed ==> $member_id"
    done
fi

# ==========================================================
# Define Experiment IDs [l_multi=true] -> esm-up2p0-gwl
# exp_ids=("esm-up2p0-gwl2p0" "esm-up2p0-gwl3p0" "esm-up2p0-gwl4p0")

# Define Experiment IDs [l_multi=true] -> esm-up2p0-gwl-dn
# exp_ids=("esm-up2p0-gwl2p0-50y-dn1p0" "esm-up2p0-gwl2p0-50y-dn2p0" "esm-up2p0-gwl3p0-50y-dn2p0" "esm-up2p0-gwl4p0-50y-dn2p0" "esm-up2p0-gwl4p0-50y-dn1p0")
# exp_ids=("esm-up2p0-gwl4p0-50y-dn1p0")


# -- Python Environment -- #
# Run this script in the env_optimesm conda virtual environment.

# if [ "$l_multi" = false ]; then
#     # -- Run NEMO Pipeline CLI -- #
#     # nemo_pipeline describe $config_file --log $log_file
#     nemo_pipeline run $config_file --log $log_file

# else
#     # Iterate over all experiment IDs:
#     for exp_id in "${exp_ids[@]}"; do
#         echo "Running ==> $exp_id"
#         # -- Updating Experiment IDs in config.toml -- #
#         sed -i "s|esm-[^/_]*|$exp_id|g" $config_file
    
#         # -- Run NEMO Pipeline CLI -- #
#         # nemo_pipeline describe $config_file --log $log_file
#         nemo_pipeline run $config_file --log $log_file
#         echo "Completed ==> $exp_id" 
#     done
# fi
