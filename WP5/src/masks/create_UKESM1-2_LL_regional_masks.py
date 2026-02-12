"""
create_UKESM1-2_LL_regional_masks.py

Description: Script to define SPNA, GIN and AO regional masks for UKESM1-2.

Date Created: 04-12-2025

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
fpath = "/g100/home/userexternal/otooth00/OptimESM/data/CINECA/MOHC/UKESM1_Ofx/domain_cfg_Ofx_UKESM1.nc"
ds_domain = xr.open_dataset(fpath).rename({"z": "nav_lev"})

# Define open-ocean mask (bathymetry > 250 m):
mask_oo = (ds_domain['bathy_metry'] > 250).rename({"y": "j", "x": "i"})

# Define path to eORCA1 monthly mean outputs:
exp_id = "r1i1p1f1"
exp_name = "esm-hist"
realm_name = "Omon"
variable_name = "thetao"
fdir = f"/g100_store/DRES_OptimESM/ESGF/prepub/mohc/20240619/CMIP6/CMIP/MOHC/UKESM1-2/{exp_name}/{exp_id}/{realm_name}"
fpaths_gridT = f"{fdir}/{variable_name}/gn/v*/{variable_name}_{realm_name}_UKESM1-2-LL_{exp_name}_{exp_id}_gn_*.nc"

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
                            "lev": "deptht",
                            "time": "time_counter"
                            })
print("Completed: Opened UKESM1-2 domain_cfg & gridT datasets.")

# Define dictionary of grid datasets defining eORCA1 parent model domain with no child/grand-child nests:
# Note: domain_cfg z-dimension is expected to be named 'nav_lev'.
datasets = {"parent":
            {"domain": ds_domain,
             "gridT": ds_gridT,
             }
            }

# Initialise a new NEMODataTree whose parent domain is zonally periodic & north-folding on F-points:
nemo = NEMODataTree.from_datasets(datasets=datasets, iperio=True, nftype="F")
print("Completed: Created UKESM1-2 NEMODataTree.")

# -- Create SPNA mask -- #
i_poly = np.array([226, 226, 230, 245, 257, 266, 277, 281, 285, 285, 285, 233, 226]) + 0.5
j_poly = np.array([275, 286, 291, 291, 289, 282, 278, 270, 262, 255, 255, 255, 275]) + 0.5

mask_SPNA = create_polygon_mask(lon_grid=nemo['gridT']['i'].expand_dims(dim={"j": nemo['gridT']['j'].size}, axis=0),
                                lat_grid=nemo['gridT']['j'].expand_dims(dim={"i": nemo['gridT']['i'].size}, axis=1),
                                lon_poly=i_poly,
                                lat_poly=j_poly,
                                dims=('j', 'i')
                                )

mask_SPNA_data = mask_SPNA.astype('int').values
boundary_SPNA = measure.find_contours(image=mask_SPNA_data, level=0.5, fully_connected="high")
print("Completed: Define Subpolar North Atlantic regional ocean mask.")

# -- Create GIN mask -- #
i_poly = np.array([bdy[1] for bdy in boundary_SPNA[0]][-55:-12])
j_poly = np.array([bdy[0] for bdy in boundary_SPNA[0]][-55:-12])
i_poly = np.concatenate([i_poly, np.array([263, 280, 294, 294, 281])]) + 0.5
# Translate SPNA boundary j by 1/2 grid spacing to prevent overlapping masks:
j_poly = np.concatenate([j_poly + 0.5, np.array([310, 310, 299, 278, 270])]) + 0.5

mask_GIN = create_polygon_mask(lon_grid=nemo['gridT']['i'].expand_dims(dim={"j": nemo['gridT']['j'].size}, axis=0),
                               lat_grid=nemo['gridT']['j'].expand_dims(dim={"i": nemo['gridT']['i'].size}, axis=1),
                               lon_poly=i_poly,
                               lat_poly=j_poly,
                               dims=('j', 'i')
                               )

mask_GIN_data = mask_GIN.astype('int').values
boundary_GIN = measure.find_contours(image=mask_GIN_data, level=0.5, fully_connected="high")
print("Completed: Define Nordic Seas regional ocean mask.")

# -- Create Arctic Ocean mask -- #
# 1. Right-Hand Side [Greenland, Davis Strait, Fram Strait etc.]
i_poly = np.concatenate([np.array([230, 230, 200, 200, 330, 330]),
                         np.array([bdy[1] for bdy in boundary_GIN[0]][-30:]),
                         np.array([bdy[1] for bdy in boundary_GIN[0]][:41]),
                         np.array([245, 230])]) + 0.5

# Translate GSR northern boundary j by 1/2 grid spacing to prevent overlapping masks:
j_poly = np.concatenate([np.array([291, 306, 306, 330, 330, 295]),
                         np.array([bdy[0] for bdy in boundary_GIN[0]][-30:]) + 0.5,
                         np.array([bdy[0] for bdy in boundary_GIN[0]][:41]) + 0.5,
                         np.array([291, 291])]) + 0.5

mask_AO_RHS = create_polygon_mask(lon_grid=nemo['gridT']['i'].expand_dims(dim={"j": nemo['gridT']['j'].size}, axis=0),
                                  lat_grid=nemo['gridT']['j'].expand_dims(dim={"i": nemo['gridT']['i'].size}, axis=1),
                                  lon_poly=i_poly,
                                  lat_poly=j_poly,
                                  dims=('j', 'i')
                                  )

mask_AO_RHS_data = mask_AO_RHS.astype('int').values
boundary_AO_RHS = measure.find_contours(image=mask_AO_RHS_data, level=0.5, fully_connected="high")

# 2. Left-Hand Side [Bering Strait etc.]
mask_AO_LHS = create_polygon_mask(lon_grid=nemo['gridT']['i'].expand_dims(dim={"j": nemo['gridT']['j'].size}, axis=0),
                                  lat_grid=nemo['gridT']['j'].expand_dims(dim={"i": nemo['gridT']['i'].size}, axis=1),
                                  lon_poly=[40, 158, 158, 40, 40],
                                  lat_poly=[285, 285, 331, 331, 285],
                                  dims=('j', 'i')
                                  )

# Combined Arctic Ocean mask:
mask_AO = (mask_AO_RHS | mask_AO_LHS)
print("Completed: Define Arctic Ocean regional ocean mask.")

# -- Write to netCDF file -- #
# Total regions:
ds_out = xr.Dataset()
ds_out["spnamsk"] = xr.DataArray(data=mask_SPNA, dims=("j", "i"))
ds_out["ginmsk"] = xr.DataArray(data=mask_GIN, dims=("j", "i"))
ds_out["aomsk"] = xr.DataArray(data=mask_AO, dims=("j", "i"))\

# Open-Ocean region component:
ds_out["spnamsk_oo"] = xr.DataArray(data=mask_SPNA & mask_oo, dims=("j", "i"))
ds_out["ginmsk_oo"] = xr.DataArray(data=mask_GIN & mask_oo, dims=("j", "i"))
ds_out["aomsk_oo"] = xr.DataArray(data=mask_AO & mask_oo, dims=("j", "i"))
# Shelf Sea region component:
ds_out["spnamsk_ss"] = xr.DataArray(data=mask_SPNA & ~mask_oo, dims=("j", "i"))
ds_out["ginmsk_ss"] = xr.DataArray(data=mask_GIN & ~mask_oo, dims=("j", "i"))
ds_out["aomsk_ss"] = xr.DataArray(data=mask_AO & ~mask_oo, dims=("j", "i"))

# Update coordinate dimensions:
ds_out = ds_out.assign_coords({"gphit": ds_domain["gphit"].rename({"y": "j", "x": "i"}),
                               "glamt": ds_domain["glamt"].rename({"y": "j", "x": "i"})
                               }
                              )
outfilepath = "/g100/home/userexternal/otooth00/OptimESM/data/CINECA/MOHC/UKESM1_Ofx/regional_masks_Ofx_UKESM1-2.nc"
ds_out.to_netcdf(outfilepath)
print(f"Completed: Saved regional [SPNA, GIN, AO] masks to netCDF file -> {outfilepath}")
