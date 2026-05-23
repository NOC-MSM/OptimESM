"""
create_UKESM1-2_LL_hemisphere_masks.py

Description: Script to define Northern & Southern Hemisphere masks for UKESM1-2-LL.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import dependencies -- #
import xarray as xr

# -- Create NEMODataTree -- #
# Define path to UKESM1-2 domain_cfg:
fpath = "/home/users/o_tooth/NOC/OptimESM/JASMIN/data/UKESM1-2-LL/Ofx/eORCA1_domain_cfg_UKESM1.nc"
ds_domain = xr.open_dataset(fpath).rename({"z": "nav_lev"})

# Define Northern Hemisphere ocean mask [latitude > 0N]:
mask_NH = (ds_domain['gphit'] > 0)
mask_NH = mask_NH.rename({"y": "j", "x": "i"})
mask_NH = mask_NH.assign_attrs({"long_name": "Northern Hemisphere ocean mask",
                                "comment": "Northern Hemisphere is defined as the region where latitude > 0N for all ocean scalar grid points in the UKESM1-2-LL domain."
                                })
# Convert to xarray.Dataset:
ds_out = mask_NH.to_dataset(name="mask_NH")

# Define Southern Hemisphere ocean mask [latitude < 0N]:
mask_SH = (ds_domain['gphit'] < 0)
mask_SH = mask_SH.rename({"y": "j", "x": "i"})
mask_SH = mask_SH.assign_attrs({"long_name": "Southern Hemisphere ocean mask",
                                "comment": "Southern Hemisphere is defined as the region where latitude < 0N for all ocean scalar grid points in the UKESM1-2-LL domain."
                                })
# Add to xarray.Dataset:
ds_out["mask_SH"] = mask_SH

# Update coordinate dimensions:
ds_out = ds_out.assign_coords({"gphit": ds_domain["gphit"].rename({"y": "j", "x": "i"}),
                               "glamt": ds_domain["glamt"].rename({"y": "j", "x": "i"})
                               }
                              )
outfilepath = "/home/users/o_tooth/NOC/OptimESM/JASMIN/data/UKESM1-2-LL/Ofx/hemisphere_masks_Ofx_UKESM1-2-LL.nc"
ds_out.to_netcdf(outfilepath)
print(f"Completed: Saved Northern & Southern Hemisphere masks to netCDF file -> {outfilepath}")
