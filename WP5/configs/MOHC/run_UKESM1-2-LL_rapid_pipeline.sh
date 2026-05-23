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
config_file=rapid/config_UKESM1_esm-hist.toml
log_file=UKESM1-2-LL_rapid_pipeline.log

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
