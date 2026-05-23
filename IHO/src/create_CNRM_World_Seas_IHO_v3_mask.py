"""
create_CNRM_World_Seas_IHO_v3_mask.py

Description:
Python script to create a 2-dimensional tracer variable labelling
T-grid points of the CNRM ocean model grid according to the IHO
sea areas they are contained within.

Created By: 
Ollie Tooth (oliver.tooth@noc.ac.uk)

"""
# -- Import required packages -- #
import logging
import numpy as np
import xarray as xr
from tqdm import tqdm
import geopandas as gpd
from matplotlib.path import Path


# -- Define utility functions -- #
def extract_mask_cells(glamt, gphit, lon_poly, lat_poly):
    """
    Extract grid cells contained within an IHO World Sea mask defined using a polygon.

    Parameters
    ----------
    glamt : DataArray
        Longitudes of model grid defined at T-points (2-dimensional array).
    gphit : DataArray
        Latitudes of model grid defined at T-points (2-dimensional array).
    lon_poly : ndarray
        Longitudes defining a polygon used to mask model grid (1-dimensional array).
    lat_poly : ndarray
        Latitudes defining a polygon used to mask model grid (1-dimensional array).

    Returns
    -------
    mask_mdl : DataArray
        Boolean mask identifying the model grid cells defined at T-points which are
        contained (True) inside the specified polyon.
    """
    # Define shape of model grid:
    mdl_shape = glamt.shape

    # Convert longitudes & latitudes to ndarrays:
    glamt = glamt.values
    gphit = gphit.values

    # Prepare polygon coordinates:
    poly_coords = list(zip(lon_poly, lat_poly))
    # Define polygon path object using coordinate tuples:
    polygon = Path(poly_coords)

    # Find coordinates of polygon bbox:
    lon_min = np.nanmin(lon_poly)
    lon_max = np.nanmax(lon_poly)
    lat_min = np.nanmin(lat_poly)
    lat_max = np.nanmax(lat_poly)
    # Determine model (i,j) coords of grid points inside bbox:
    ind_bbox = np.argwhere((glamt >= lon_min) & (glamt <= lon_max) & (gphit >= lat_min) & (gphit <= lat_max))
    i_bbox = ind_bbox[:, 1]
    j_bbox = ind_bbox[:, 0]

    if i_bbox.size > 0:
        # Prepare NEMO model coordinates in bbox:
        lon_mdl = glamt[j_bbox, i_bbox]
        lat_mdl = gphit[j_bbox, i_bbox]
        mdl_coords = np.array(list(zip(lon_mdl, lat_mdl)))

        # Create boolean mask for enclosed NEMO model grid cells and reshape to 2d:
        mask_points = polygon.contains_points(mdl_coords)
        # Define initial model grid cell mask ndarray [False]:
        mask_mdl = np.zeros(mdl_shape, dtype=np.bool)
        # Identify where model grid cells are inside polygon
        mask_mdl[j_bbox[mask_points], i_bbox[mask_points]] = True

    else: 
        # Return model grid cell mask ndarray [False]:
        mask_mdl = np.zeros(mdl_shape, dtype=np.bool)

    return xr.DataArray(data=mask_mdl, dims=['y', 'x'])


def main(
    domain_fpath: str,
    shp_fpath: str,
    out_fpath: str
):
    # -- Open NEMO model Domain and Bathymetry Files -- #
    ds_domain = xr.open_dataset(domain_fpath)[['longitude', 'latitude']]

    # -- Open IHO sea masks shapefile as GeoDataFrame -- #
    df_iho = gpd.read_file(shp_fpath)
    df_iho['IHO'] = df_iho.index

    # -- Define NEMO model inputs -- #
    glamt = ds_domain['longitude'].squeeze()
    gphit = ds_domain['latitude'].squeeze()
    logging.info('Completed: Prepared Model Grid Data')

    # -- Define empty tracer array to store IHO sea mask IDs -- #
    iho_tracer = np.zeros(glamt.shape)
    iho_tracer[:, :] = np.nan
    iho_tracer = xr.DataArray(data=iho_tracer, dims=['y', 'x'])
    logging.info('Completed: Created IHO 2-dimensional Tracer Variable.')

    # -- Find NEMO model grid cells contained in each IHO sea mask -- #
    # Define no. IHO entries:
    n_entries = df_iho.index.size

    # Define empty list to store IHO IDs & Names:
    iho_ids = []
    iho_MRGIDs = []
    iho_names = []

    # Iterate over World Seas IHO v3 entries:
    logging.info('In Progress: Extracting IHO World Sea Areas masks from NEMO domain.')
    for iho_n in tqdm(range(n_entries)):
        # Define IHO name:
        iho_name = df_iho.loc[df_iho.index == iho_n]['NAME'].item()
        # Store IHO ID and Name:
        iho_ids.append(iho_n)
        iho_names.append(iho_name)
        # Store IHO MRGID:
        iho_MRGIDs.append(df_iho.loc[df_iho.index == iho_n]['MRGID'].item())

        # -- 1. Extract IHO World Sea Area with MultiPolygon coordinates -- #
        if 'MultiPolygon' in df_iho.loc[df_iho.index == iho_n]['geometry'].geom_type.values:
            # Explode MultiPolygon to DataFrame of Polygons xy coordinates:
            df_poly = df_iho.loc[df_iho.index == iho_n].explode(index_parts=True)['geometry'].get_coordinates()
            # Get the number of polygons:
            n_poly = df_poly.index.to_series().nunique()

            # Iterate over polygons:
            for poly_n in range(n_poly):
                # Extract closed polygon coordinates:
                lon_poly = df_poly.loc[(iho_n, poly_n)]['x'].values
                lat_poly = df_poly.loc[(iho_n, poly_n)]['y'].values

                # Determine NEMO model grid cell mask:
                mask_traj = extract_mask_cells(glamt, gphit, lon_poly, lat_poly)
                # Ignore masks with no grid cells in IHO:
                if mask_traj.sum() > 0:
                    # Update IHO IDs stored in 2d tracer variable:
                    iho_tracer = xr.where(mask_traj, np.int32(iho_n), iho_tracer)

        # -- 2. Extract IHO World Sea Area with Polygon coordinates -- #
        else:
            # Collect DataFrame of Polygon xy coordinates:
            df_poly = df_iho.loc[df_iho.index == iho_n]['geometry'].get_coordinates()
            # Extract polygon coordinates:
            lon_poly = df_poly['x'].values
            lat_poly = df_poly['y'].values

            # Select only the (lon, lat) coordinates required to close the polygon:
            lon_inds = np.argwhere(lon_poly == lon_poly[0])
            lat_inds = np.argwhere(lat_poly == lat_poly[0])
            # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
            n_closed = [n for n in lon_inds if n in lat_inds][1].item()
            if n_closed > 1:
                # Select only the coordinates required to close the polygon:
                lon_poly = lon_poly[:n_closed]
                lat_poly = lat_poly[:n_closed]
            else:
                raise ValueError("IHO World Sea Area polygon is not closed.")

            # Determine NEMO model grid cell mask:
            mask_traj = extract_mask_cells(glamt, gphit, lon_poly, lat_poly)
            # Ignore masks with no grid cells in IHO:
            if mask_traj.sum() > 0:
                # Update IHO IDs stored in 2d tracer variable:
                iho_tracer = xr.where(mask_traj, np.int32(iho_n), iho_tracer)

    logging.info('Completed: Identified IHO World Sea Area mask on NEMO model grid.')

    # -- Saving IHO Tracer DataArray to File -- #
    # Add IHO IDs, Names and MRGIDs as DataArrays:
    ds_iho = iho_tracer.to_dataset(name='mask_iho')
    ds_iho['name_iho'] = xr.DataArray(data=iho_names, dims=['id'], coords={'id': iho_ids})
    ds_iho['mrgid_iho'] = xr.DataArray(data=iho_MRGIDs, dims=['id'], coords={'id': iho_ids})
    ds_iho.to_netcdf(out_fpath, unlimited_dims="time")

    logging.info('Completed: Saved IHO World Sea Area masks on NEMO model grid to netCDF file.')

if __name__ == "__main__":

    # --- Configure Logging --- #
    logging.basicConfig(
        filename="CNRM-ESM2-1_World_Seas_IHO_masks.log",
        encoding="utf-8",
        filemode="a",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
        level=logging.INFO,
        )

    # --- Run Main Function --- #
    main(
        domain_fpath="/g100/home/userexternal/otooth00/OptimESM/data/CINECA/CNRM/CNRM_Ofx/mesh_mask_Ofx_CNRM.nc",
        shp_fpath="/g100_work/optim_IAC/research/noc/otooth/OptimESM/IHO/World_Seas_IHO_v3/World_Seas_IHO_v3_combined.shp",
        out_fpath="/g100_work/optim_IAC/research/noc/otooth/OptimESM/IHO/masks/CNRM-ESM2-1_World_Seas_IHO_v3_mask.nc"
    )
