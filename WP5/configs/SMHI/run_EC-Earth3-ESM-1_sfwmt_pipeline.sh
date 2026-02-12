#!/bin/bash

set -euo pipefail
# ================================================================
# run_EC-Earth3-ESM-1_pipeline.sh
#
# Description: Run EC-Earth3-ESM-1 NEMO Pipeline in current process.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-11-11   
# ================================================================

# -- Input arguments to NEMO Pipeline -- #
# Define filepaths:
config_file=siarea/config_EC-Earth3-ESM-1_esm-up2p0-gwl-dn.toml
log_file=EC-Earth3-ESM-1_siarea_pipeline.log

# -- Python Environment -- #
# Run this script in the env_optimesm conda virtual environment.

# -- Run NEMO Pipeline CLI -- #
# nemo_pipeline describe $config_file --log $log_file
nemo_pipeline run $config_file --log $log_file
