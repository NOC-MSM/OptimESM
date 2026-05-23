"""
create_CNRM-ESM2-1_domain_cfg.py

Description: Script to define complete NEMO Cookbook compatible
domain_cfg.nc file for CNRM-ESM2-1.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import dependencies -- #
import numpy as np
import xarray as xr

# -- Create domain_cfg.nc -- #
# Define path to domain_cfg:
fpath = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/CNRM-ESM2-1/Ofx/eORCA1L75_CNRM-ESM2-1.nc"
ds_domain = xr.open_dataset(fpath).rename({"z": "nav_lev"}).drop_vars(["AREA", "VOLUME", "nav_lev"])

# Define path to gridW file (wmask & wmaskutil definition):
fpath = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/CNRM-ESM2-1/esm-up2p0/r1i1p2f2/Omon/wo/gn/wo_Omon_CNRM-ESM2-1_02Kpd_r1i1p2f2_gn_185001-187412.nc"
ds_gridW = xr.open_dataset(fpath)

# wmask from wo:
wo = ds_gridW['wo'][0, :, :, :].squeeze(drop=True).drop_vars("time").rename({"lev": "nav_lev"})
wmask = 2 + ~np.isnan(wo).drop_vars(["nav_lev", "lat", "lon"]).astype('int8')
wmask.name = "wmask"
wmask.attrs = {}

# wmaskutil from wmask:
wmaskutil = wmask[0, :, :]
wmaskutil.name = "wmaskutil"
wmaskutil.attrs = {}

# Assign new variables to domain_cfg:
ds_domain["wmask"] = wmask
ds_domain["wmaskutil"] = wmaskutil

# Write to local netCDF file:
outfpath = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/CNRM-ESM2-1/Ofx/eORCA1L75_domain_cfg_CNRM-ESM2-1.nc"
ds_domain.to_netcdf(outfpath)
print(f"Completed: Saved CNRM-ESM2-1 domain_cfg to netCDF file -> {outfpath}")
