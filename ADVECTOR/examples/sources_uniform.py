import time 
from pathlib import Path as path

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union
from global_land_mask import globe

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

"""
Package for land mask is: https://github.com/toddkarin/global-land-mask
Coordinates of North Pacific Ocean https://www.marineregions.org/gazetteer.php/gazetteer.php?p=details&id=1908 
"""

check = False

start = time.time()
wesn = list([117, -78, 0, 67]) 
wesn = [wesn[i] if wesn[i] >= 0 else wesn[i]+360 for i in range(4)]
res_km = 10 # km

R = 6371 #km
deg_km = 6371*np.pi/180 # km per degree
dlat = res_km/deg_km 
lats = np.arange(wesn[2], wesn[3], dlat)
dlon = res_km / (deg_km * np.cos(np.deg2rad(lats)))

lons = []
for step in dlon:
    l = np.arange(wesn[0], wesn[1], step)
    l[l > 180] -= 360
    lons.append(l)
# lons = [np.arange(wesn[0], wesn[1], step) for step in dlon]

max_len = max(len(lon) for lon in lons)
lon_grid = np.full((len(lats), max_len), np.nan)
lat_grid = np.full((len(lats), max_len), np.nan)  
for i, lon in enumerate(lons):
    lon_grid[i, :len(lon)] = lon  
    lat_grid[i, :len(lon)] = lats[i]
end = time.time()
print(f"Time taken to create grid: {end-start}")

start = time.time()
globe_land_mask = globe.is_land(lat_grid, lon_grid)
end = time.time()
print(f"Time taken to create land mask: {end-start}")

sources_df = xr.open_dataset(path("/Volumes/storage4/ADVECTOR/sources/uniform/T1-500_RD2020.nc"))  
lon_grid_masked = lon_grid[~globe_land_mask]
lat_grid_masked = lat_grid[~globe_land_mask]
release_date = np.ones(lon_grid_masked.shape[0])
release_date = np.repeat(pd.to_datetime("2020-01-01 00:00:00", format="%Y-%m-%d %H:%M:%S"), lon_grid_masked.shape[0])
sources_xr = xr.Dataset()
sources_xr["lon"] = xr.DataArray(lon_grid_masked, dims=("p_id",), coords={"p_id": np.arange(lon_grid_masked.shape[0])})
sources_xr["lat"] = xr.DataArray(lat_grid_masked, dims=("p_id",), coords={"p_id": np.arange(lat_grid_masked.shape[0])})
sources_xr["release_date"] = xr.DataArray(release_date, dims=("p_id",), coords={"p_id": np.arange(lon_grid_masked.shape[0])})
# sources_xr['unsd'] = xr.DataArray(np.ones(lon_grid_masked.shape[0]), dims=("p_id",), coords={"p_id": np.arange(lon_grid_masked.shape[0])})
# add attributes
sources_xr.attrs["title"] = "Particle sources for ADVECTOR-type Lagrangian Dispersion Model"
sources_xr.attrs["institution"] = "The Ocean Cleanup"
sources_xr.attrs["source"] = "Uniform grid with 10 km resolution over the North Pacific Ocean"
sources_xr.to_netcdf(path("/Volumes/storage4/ADVECTOR/sources/uniform/T0-1e6_RD2020_land.nc"))
# plot 
fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': ccrs.PlateCarree()})
ax.coastlines()
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
sources_xr.plot.scatter(x="lon", y="lat", ax=ax, color='r', s=1)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
plt.savefig(path("/Volumes/storage4/ADVECTOR/sources/uniform/T0-1e6_RD2020_land.png"))

# url = "https://naciscdn.org/naturalearth/110m/physical/ne_110m_coastline.zip"
# url = "https://naciscdn.org/naturalearth/110m/physical/ne_110m_land.zip"
# land = gpd.read_file(url)
# land.boundary.plot(ax=ax)

if check == True:
    start = time.time()
    points = [Point(lon, lat) for lon, lat in zip(lon_grid.flatten(), lat_grid.flatten())]
    points_gdf = gpd.GeoDataFrame(geometry=points, crs="EPSG:4326")
    points_within_land = [Point(lon, lat) for lon, lat in zip(lon_grid_masked.flatten(), lat_grid_masked.flatten())]
    points_within_land_gdf = gpd.GeoDataFrame(geometry=points_within_land, crs="EPSG:4326")
    # masked = points_gdf.within(land.geometry.unary_union)
    end = time.time()
    print(f"Time taken to create GeoDataFrame: {end-start}")

    fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    points_gdf.plot(ax=ax, color='r', markersize=1)
    points_within_land_gdf.plot(ax=ax, color='b', markersize=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
    plt.show()

