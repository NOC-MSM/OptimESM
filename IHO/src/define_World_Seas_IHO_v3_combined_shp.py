"""
define_World_Seas_IHO_v3_combined_shp.py

Description:
Python script to define an updated shapefile containing MultiPolygon geometries
IHO sea areas crossing the international dateline, including the Bering Sea,
North Pacific Ocean, and South Pacific Ocean in the IHO sea areas.

Created By: 
Ollie Tooth (oliver.tooth@noc.ac.uk)

"""
# -- Import required packages -- #
import logging
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon


def main(
    shp_fpath: str,
    out_fpath: str
    ):
    logging.info("In Progress: Defining World Seas IHO v3 combined shapefile")

    # -- Open IHO sea masks shapefile as GeoDataFrame -- #
    df_iho = gpd.read_file(shp_fpath)
    logging.info(f"Completed: Opened World Seas IHO v3 shapefile from: {shp_fpath}")

    # ================== North Pacific Ocean ================== #
    df_poly = df_iho.loc[df_iho.NAME == 'North Pacific Ocean']['geometry'].get_coordinates()
    lon_poly = df_poly['x'].values
    lat_poly = df_poly['y'].values

    # Polygon 1 - West:
    n_start = 0 # start index of western polygon
    lon_inds = np.argwhere(lon_poly == lon_poly[n_start])
    lat_inds = np.argwhere(lat_poly == lat_poly[n_start])
    # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
    n_closed = [n for n in lon_inds if n in lat_inds][1].item()
    if n_closed > 1:
        # Select only the coordinates required to close the polygon:
        p1 = Polygon(list(zip(lon_poly[n_start:n_closed], lat_poly[n_start:n_closed])))

    # Polygon 2 - East:
    n_start = 310830 # starting index for the eastern polygon
    lon_inds = np.argwhere(lon_poly == lon_poly[n_start])
    lat_inds = np.argwhere(lat_poly == lat_poly[n_start])
    # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
    n_closed = [n for n in lon_inds if n in lat_inds][1].item()
    if n_closed > 1:
        # Select only the coordinates required to close the polygon:
        p2 = Polygon(list(zip(lon_poly[n_start:n_closed], lat_poly[n_start:n_closed])))

    # Replace existing Polygon Geometry with MultiPolygon:
    df_iho.loc[df_iho.NAME == 'North Pacific Ocean', "geometry"] = MultiPolygon([p1, p2])

    logging.info("Completed: Defined MultiPolygon geometry for North Pacific Ocean")

    # ================== South Pacific Ocean ================== #
    df_poly = df_iho.loc[df_iho.NAME == 'South Pacific Ocean']['geometry'].get_coordinates()
    lon_poly = df_poly['x'].values
    lat_poly = df_poly['y'].values

    # Polygon 1 - West:
    n_start = 0 # start index of western polygon
    lon_inds = np.argwhere(lon_poly == lon_poly[n_start])
    lat_inds = np.argwhere(lat_poly == lat_poly[n_start])
    # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
    n_closed = [n for n in lon_inds if n in lat_inds][1].item()
    if n_closed > 1:
        # Select only the coordinates required to close the polygon:
        p1 = Polygon(list(zip(lon_poly[n_start:n_closed], lat_poly[n_start:n_closed])))

    # Polygon 2 - East:
    n_start = 280157 # starting index for the eastern polygon
    lon_inds = np.argwhere(lon_poly == lon_poly[n_start])
    lat_inds = np.argwhere(lat_poly == lat_poly[n_start])
    # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
    n_closed = [n for n in lon_inds if n in lat_inds][1].item()
    if n_closed > 1:
        # Select only the coordinates required to close the polygon:
        p2 = Polygon(list(zip(lon_poly[n_start:n_closed], lat_poly[n_start:n_closed])))

    # Replace existing Polygon Geometry with MultiPolygon:
    df_iho.loc[df_iho.NAME == 'South Pacific Ocean', "geometry"] = MultiPolygon([p1, p2])

    logging.info("Completed: Defined MultiPolygon geometry for South Pacific Ocean")

    # ================== Bering Sea ================== #
    df_poly = df_iho.loc[df_iho.NAME == 'Bering Sea']['geometry'].get_coordinates()
    lon_poly = df_poly['x'].values
    lat_poly = df_poly['y'].values

    # Polygon 1 - West:
    n_start = 0 # start index of western polygon
    lon_inds = np.argwhere(lon_poly == lon_poly[n_start])
    lat_inds = np.argwhere(lat_poly == lat_poly[n_start])
    # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
    n_closed = [n for n in lon_inds if n in lat_inds][1].item()
    if n_closed > 1:
        # Select only the coordinates required to close the polygon:
        p1 = Polygon(list(zip(lon_poly[n_start:n_closed], lat_poly[n_start:n_closed])))

    # Polygon 2 - East:
    n_start = 20800 # starting index for the eastern polygon
    lon_inds = np.argwhere(lon_poly == lon_poly[n_start])
    lat_inds = np.argwhere(lat_poly == lat_poly[n_start])
    # Use 2nd index where (lon, lat) == (lon_ini, lat_ini) to subset closed polygon:
    n_closed = [n for n in lon_inds if n in lat_inds][1].item()
    if n_closed > 1:
        # Select only the coordinates required to close the polygon:
        p2 = Polygon(list(zip(lon_poly[n_start:n_closed], lat_poly[n_start:n_closed])))

    # Replace existing Polygon Geometry with MultiPolygon:
    df_iho.loc[df_iho.NAME == 'Bering Sea', "geometry"] = MultiPolygon([p1, p2])

    logging.info("Completed: Defined MultiPolygon geometry for Bering Sea")

    # -- Write updated World Seas IHO v3 GeoDataFrame to shapefile -- #
    df_iho.to_file(out_fpath, driver='ESRI Shapefile')

    logging.info(f"Completed: Saved Updated World Seas IHO v3 shapefile to: {out_fpath}")

if __name__ == "__main__":
    # --- Configure Logging --- #
    logging.basicConfig(
        filename="World_Seas_IHO_v3_combined.log",
        encoding="utf-8",
        filemode="a",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
        level=logging.INFO,
        )

    # --- Run Main Function --- #
    main(
        shp_fpath="/g100_work/optim_IAC/research/noc/otooth/OptimESM/IHO/World_Seas_IHO_v3/World_Seas_IHO_v3.shp",
        out_fpath="/g100_work/optim_IAC/research/noc/otooth/OptimESM/IHO/World_Seas_IHO_v3/World_Seas_IHO_v3_combined.shp"
    )
