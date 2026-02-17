"""
create_IPSL-CM6-ESMCO2_SO_regional_masks.py

Description: Script to define Southern Ocean regional masks for IPSL-CM6-ESMCO2.

Date Created: 15-02-2026

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import dependencies -- #
import numpy as np
import xarray as xr
from skimage import measure
from nemo_cookbook import NEMODataTree
from nemo_cookbook.masks import create_polygon_mask

# -- Create NEMODataTree -- #
# Define path to UKESM1-2 domain_cfg:
fpath = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/IPSL-CM6-ESMCO2/Ofx/eORCA1.2_mesh_mask_IPSL.nc"
ds_domain = xr.open_dataset(fpath).rename({"z": "nav_lev"}).squeeze()

fpath_bathy = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/IPSL-CM6-ESMCO2/Ofx/eORCA_R1_bathy_meter_v2.2_IPSL.nc"
ds_bathy = xr.open_dataset(fpath_bathy).squeeze()

# Define open-ocean mask (bathymetry > 250 m):
mask_oo = (ds_bathy['Bathymetry'] > 250).rename({"y": "j", "x": "i"})

# Define path to eORCA1 monthly mean outputs:
exp_id="r1i1p3f1"
exp_name="esm-hist"
realm_name="Omon"
variable_name="thetao"
fdir = f"/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/IPSL-CM6-ESMCO2/{exp_name}/{exp_id}/{realm_name}"
fpaths_gridT=f"{fdir}/{variable_name}/gn/v*/{variable_name}_{realm_name}_IPSL-CM6-ESMCO2_{exp_name}_{exp_id}_gn_*.nc"

# Define CFDatetimeCoder to decode time coords:
coder = xr.coders.CFDatetimeCoder(time_unit="s")
ds_gridT = xr.open_mfdataset(fpaths_gridT,
                             data_vars="minimal",
                             compat="no_conflicts",
                             decode_times=coder,
                             parallel=False,
                             engine="netcdf4"
                            )

ds_gridT = ds_gridT.rename({"i": "x",
                            "j": "y",
                            "lev": "deptht",
                             "time": "time_counter"
                              })
print("Completed: Opened IPSL-CM6-ESMCO2 domain_cfg & gridT datasets.")

# Define dictionary of grid datasets defining eORCA1 parent model domain with no child/grand-child nests:
# Note: domain_cfg z-dimension is expected to be named 'nav_lev'.
datasets = {"parent":
            {"domain": ds_domain,
             "gridT": ds_gridT,
             }
            }

# Initialise a new NEMODataTree whose parent domain is zonally periodic & north-folding on F-points:
nemo = NEMODataTree.from_datasets(datasets=datasets, iperio=True, nftype="F", read_mask=True)
print("Completed: Created IPSL-CM6-ESMCO2 NEMODataTree.")

# -- Create Southern Ocean sector masks -- #
# Weddell Sea:
lon_WS = [-65, 20, 20, -65, -65]
lat_WS = [-50, -50, -90, -90, -50]
mask_WS = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_WS, lat_poly=lat_WS)

# Indian Ocean:
lon_IO = [20, 90, 90, 20, 20]
lat_IO = [-50, -50, -90, -90, -50]
mask_IO = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_IO, lat_poly=lat_IO)

# Western Pacific:
lon_WP = [90, 160, 160, 90, 90]
lat_WP = [-50, -50, -90, -90, -50]
mask_WP = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_WP, lat_poly=lat_WP)

# Ross Sea:
lon_RS_LHS = [-180, -130, -130, -180, -180]
lat_RS_LHS = [-50, -50, -90, -90, -50]
lon_RS_RHS = [160, 180, 180, 160, 160]
lat_RS_RHS = [-50, -50, -90, -90, -50]
mask_RS_LHS = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_RS_LHS, lat_poly=lat_RS_LHS)
mask_RS_RHS = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_RS_RHS, lat_poly=lat_RS_RHS)
mask_RS = (mask_RS_LHS | mask_RS_RHS)

# Amundsen-Bellingshausen Sea:
lon_ABS = [-130, -65, -65, -130, -130]
lat_ABS = [-50, -50, -90, -90, -50]
mask_ABS = nemo.mask_with_polygon(grid="gridT", lon_poly=lon_ABS, lat_poly=lat_ABS)

print("Completed: Define Southern Ocean regional ocean masks.")

# -- Write to netCDF file -- #
# Total regions:
ds_out = xr.Dataset()
ds_out["wsmsk"] = xr.DataArray(data=mask_WS, dims=("j", "i"))
ds_out["iomsk"] = xr.DataArray(data=mask_IO, dims=("j", "i"))
ds_out["wpmsk"] = xr.DataArray(data=mask_WP, dims=("j", "i"))
ds_out["rsmsk"] = xr.DataArray(data=mask_RS, dims=("j", "i"))
ds_out["absmsk"] = xr.DataArray(data=mask_ABS, dims=("j", "i"))

# Open-Ocean region component:
ds_out["wsmsk_oo"] = xr.DataArray(data=mask_WS & mask_oo, dims=("j", "i"))
ds_out["iomsk_oo"] = xr.DataArray(data=mask_IO & mask_oo, dims=("j", "i"))
ds_out["wpmsk_oo"] = xr.DataArray(data=mask_WP & mask_oo, dims=("j", "i"))
ds_out["rsmsk_oo"] = xr.DataArray(data=mask_RS & mask_oo, dims=("j", "i"))
ds_out["absmsk_oo"] = xr.DataArray(data=mask_ABS & mask_oo, dims=("j", "i"))

# Shelf Sea region component:
ds_out["wsmsk_ss"] = xr.DataArray(data=mask_WS & ~mask_oo, dims=("j", "i"))
ds_out["iomsk_ss"] = xr.DataArray(data=mask_IO & ~mask_oo, dims=("j", "i"))
ds_out["wpmsk_ss"] = xr.DataArray(data=mask_WP & ~mask_oo, dims=("j", "i"))
ds_out["rsmsk_ss"] = xr.DataArray(data=mask_RS & ~mask_oo, dims=("j", "i"))
ds_out["absmsk_ss"] = xr.DataArray(data=mask_ABS & ~mask_oo, dims=("j", "i"))

# Update coordinate dimensions:
ds_out = ds_out.assign_coords({"gphit": ds_domain["gphit"].rename({"y": "j", "x": "i"}),
                               "glamt": ds_domain["glamt"].rename({"y": "j", "x": "i"})
                               }
                              )
outfilepath = outfilepath = "/g100_work/optim_IAC/research/noc/otooth/OptimESM/data/IPSL-CM6-ESMCO2/Ofx/SO_regional_masks_Ofx_IPSL-CM6-ESMCO2.nc"
ds_out.to_netcdf(outfilepath)
print(f"Completed: Saved Southern Ocean regional masks to netCDF file -> {outfilepath}")
