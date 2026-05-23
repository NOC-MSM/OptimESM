# Preparing World Seas IHO version 3 Geometries for Masking OGCMs

### Step 1: Download IHO Sea Areas form marineregions.org

https://www.marineregions.org/downloads.php

### Step 2: Define combined World Seas Areas IHO v3

- World Seas Areas v3 which cross the international dateline need to be transformed from a single intersecting Polygon (non-Esri compliant) to a MultiPolygon geometry.

- To do this, we need to run the `define_World_Seas_IHO_v3_combined_shp.py` script as follows:

```python
python3 define_World_Seas_IHO_v3_combined_shp.py
```

### Step 3: Create World Seas Areas IHO v3 masks for given Ocean General Circulation Model

- Create World Seas Area IHO v3 masks for our given Ocean General Circulation Model. 
- This produces a netCDF file containing three DataArrays **mask_iho** - integer IDs identifying the World Sea Area to which each model grid cell belongs, **name_iho** - name of World Sea Area corresponding to each integer ID,  **mrgid_iho** - MRGID value corresponding to World Sea Area.
- To do this, we need to run the `create_World_Seas_IHO_v3_masks.py` script as follows:

```python
python3 create_World_Seas_IHO_v3_masks.py
```