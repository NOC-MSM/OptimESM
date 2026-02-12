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
config_file=si/config_UKESM1_esm-up2p0-gwl-dn.toml
log_file=UKESM1-2-LL_si_pipeline.log

# -- Python Environment -- #
# Run this script in the env_optimesm conda virtual environment.

# -- Run NEMO Pipeline CLI -- #
# nemo_pipeline describe $config_file --log $log_file
nemo_pipeline run $config_file --log $log_file
