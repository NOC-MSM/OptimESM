#!/bin/bash

set -euo pipefail
# ================================================================
# run_EC-Earth3-ESM-1_pipeline.sh
#
# Description: Run EC-Earth3-ESM-1 NEMO Pipeline in current process.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-05-08
# ================================================================

# -- Input arguments to NEMO Pipeline -- #
# Define filepaths:
config_file=fram/config_EC-Earth3-ESM-1_esm-hist.toml
log_file=EC-Earth3-ESM-1_fram_pipeline.log

# Run multiple pipelines:
l_multi=true

# Define Experiment IDs [l_multi=true] -> esm-hist ensemble members...
# member_ids=("r1i1p1f1" "r2i1p1f1" "r3i1p1f1" "r5i1p1f1")
member_ids=("r1i1p1f1" "r2i1p1f1")

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

