"""
create_UKESM1-2_LL_SPNA_ML_masks.py

Description: Script to define SPNA mixed layer masks for UKESM1-2-LL.

Date Created: 16-03-2026

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import dependencies -- #
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from nemo_cookbook import NEMODataTree

# -- Create NEMODataTree -- #
# Define path to UKESM1-2 domain_cfg:
fpath = "/g100/home/userexternal/otooth00/OptimESM/data/CINECA/MOHC/UKESM1_Ofx/domain_cfg_Ofx_UKESM1.nc"
ds_domain = xr.open_dataset(fpath).rename({"z": "nav_lev"})

# Define path to eORCA1 monthly mean outputs:
exp_id = "r1i1p1f1"
exp_name = "esm-piControl"
realm_name = "Omon"
variable_name = "mlotst"
fdir = f"/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/UKESM1-2-LL/{exp_name}/{exp_id}/{realm_name}"
fpaths_gridT = f"{fdir}/{variable_name}/gn/v20241002/{variable_name}_{realm_name}_UKESM1-2-LL_{exp_name}_{exp_id}_gn_*.nc"

# Define CFDatetimeCoder to decode time coords:
coder = xr.coders.CFDatetimeCoder(time_unit="s")
ds_gridT = xr.open_mfdataset(fpaths_gridT,
                             data_vars="minimal",
                             compat="no_conflicts",
                             decode_times=coder,
                             parallel=False,
                             engine="netcdf4"
                             )
# Update coord dims to standard NEMO dims:
ds_gridT = ds_gridT.rename({"i": "x",
                            "j": "y",
                            "time": "time_counter"
                            })
print("Completed: Opened UKESM1-2-LL domain_cfg & gridT datasets.")

# Define dictionary of grid datasets defining eORCA1 parent model domain with no child/grand-child nests:
# Note: domain_cfg z-dimension is expected to be named 'nav_lev'.
datasets = {"parent": {"domain": ds_domain, "gridT": ds_gridT}}

# Initialise a new NEMODataTree whose parent domain is zonally periodic & north-folding on F-points:
nemo = NEMODataTree.from_datasets(datasets=datasets, iperio=True, nftype="F")
print("Completed: Created UKESM1-2-LL NEMODataTree.")

# -- Create Labrador Sea MLD mask -- #
# Compute March mean mixed layer depth for UKESM1-2-LL esm-piControl:
mlotst_mmean = nemo['gridT/mlotst'].sel(time_counter=nemo['gridT']['time_counter'].dt.month.isin([3])).mean(dim='time_counter')
cs = plt.contour(ds_gridT['longitude'], ds_gridT['latitude'], mlotst_mmean, levels=[500])

# Extract the closed 500-m MLD polygon for the Labrador Sea:
ml_poly = cs.allsegs[0][3]
lon_LabSea_poly = np.array([coord[0] for coord in ml_poly], dtype=np.float32)
lat_LabSea_poly = np.array([coord[1] for coord in ml_poly], dtype=np.float32)
# Define mask:
mask_LabSea = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_LabSea_poly, lat_poly=lat_LabSea_poly)
mask_LabSea = mask_LabSea.astype('int').values
print("Completed: Defined Labrador Sea mixed layer [500 m MLD contour] ocean mask.")

# -- Create Irminger Sea MLD mask -- #
# Extract the closed 500-m MLD polygon for the Irminger Sea:
ml_poly = cs.allsegs[0][4]
lon_IrmSea_poly = np.array([coord[0] for coord in ml_poly], dtype=np.float32)
lat_LabSea_poly = np.array([coord[1] for coord in ml_poly], dtype=np.float32)
# Define mask:
mask_IrmSea = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_IrmSea_poly, lat_poly=lat_LabSea_poly)
mask_IrmSea = mask_IrmSea.astype('int').values
print("Completed: Defined Irminger Sea mixed layer [500 m MLD contour] ocean mask.")

# -- Write to netCDF file -- #
# Total regions:
ds_out = xr.Dataset()
ds_out["irmseamsk"] = xr.DataArray(data=mask_IrmSea, dims=("j", "i"))
ds_out["labseamsk"] = xr.DataArray(data=mask_LabSea, dims=("j", "i"))

# Update coordinate dimensions:
ds_out = ds_out.assign_coords({"gphit": ds_domain["gphit"].rename({"y": "j", "x": "i"}),
                               "glamt": ds_domain["glamt"].rename({"y": "j", "x": "i"})
                               }
                              )
outfilepath = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/UKESM1-2-LL/Ofx/MLD_masks_Ofx_UKESM1-2-LL.nc"
ds_out.to_netcdf(outfilepath)
print(f"Completed: Saved regional [Labrador Sea, Irminger Sea] MLD masks to netCDF file -> {outfilepath}")
