import tqdm     
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

def land_masker(wesn, n_cells_deg, buffer):
    """ 
    DOES NOT WORK 
    """
    
    wesn_land = [wesn[0]-10, wesn[1]+10, wesn[2]-10, wesn[3]+10]

    # Load the Natural Earth dataset for land boundaries
    url = "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip"
    land = gpd.read_file(url)
    land = land[land['geometry'].notnull()]
    land = land.cx[wesn_land[0]:wesn_land[1], wesn_land[2]:wesn_land[3]]

    lats = np.linspace(wesn[2], wesn[3], (wesn[3]-wesn[2])*n_cells_deg+1)
    lons = np.linspace(wesn[0], wesn[1], (wesn[1]-wesn[0])*n_cells_deg+1)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    points = [Point(x, y) for x, y in zip(lon_grid.flatten(), lat_grid.flatten())]
    points_gdf = gpd.GeoDataFrame(geometry=points, crs="EPSG:4326")
    points_gdf['lon'] = lon_grid.flatten()
    points_gdf['lat'] = lat_grid.flatten()

    land_buffered = land[['geometry', 'LABELRANK']].copy()
    land_buffered['geometry'] = land_buffered.buffer(buffer)

    points_gdf = gpd.sjoin(points_gdf, land_buffered, how='left', predicate='within')
    mask_gdf = points_gdf.dropna()[['lon', 'lat']]

    mask_ds = xr.Dataset(coords={'lat': lats, 'lon': lons})
    mask_ds['mask'] = xr.DataArray(np.full((len(lats), len(lons)), True), coords=[lats, lons], dims=['lat', 'lon'])
    for i in range(len(mask_gdf)):
        lon = mask_gdf.lon.values[i]
        lat = mask_gdf.lat.values[i]
        mask_ds['mask'].loc[lat, lon] = False

    return mask_ds

def advector_gridding(advector_ds, wesn, n_cells_deg):  

    res = 1/(2*n_cells_deg) 

    lats = np.linspace(wesn[2], wesn[3], (wesn[3]-wesn[2])*n_cells_deg+1)
    lons = np.linspace(wesn[0], wesn[1], (wesn[1]-wesn[0])*n_cells_deg+1)

    lats_edges = np.linspace(wesn[2]-res, wesn[3]+res, (wesn[3]-wesn[2])*n_cells_deg+2)
    lons_edges = np.linspace(wesn[0]-res, wesn[1]+res, (wesn[1]-wesn[0])*n_cells_deg+2)

    # compute number of objects in each grid cell 
    ds = xr.Dataset(coords={'time': advector_ds.time, 'lat': lats, 'lon': lons})
    ds['stats'] = xr.DataArray(np.zeros((len(advector_ds.time), len(lats), len(lons))), coords=[advector_ds.time, lats, lons], dims=['time', 'lat', 'lon'])

    # QC hist
    """ lats = np.arange(0, 2, 1)
    lons = np.arange(0, 2, 1)
    lat_bins = np.arange(-0.5, 2.5, 1) 
    lon_bins = np.arange(-0.5, 2.5, 1)
    values = [[0, 1], [1, 0]]
    mask_ds = xr.Dataset(coords={'lat': lats, 'lon': lons}, data_vars={'mask': (('lat', 'lon'), values)})
    hist, xedges, yedges = np.histogram2d(values[0], values[1], bins=(lon_bins, lat_bins))
    hist2, xedges2, yedges2 = np.histogram2d(values[0], values[1], bins=(lons, lats))  
    mask_ds['mask2'] = xr.DataArray(hist, coords=[lats, lons], dims=['lat', 'lon'])
    print(mask_ds.mask)
    print(mask_ds.mask2)
    print(hist2) 
    """
    
    # for i in tqdm(range(len(advector_ds.time))):
    for i in range(len(advector_ds.time)):
        ds_t = advector_ds.sel(time=advector_ds.time[i])

        lons = ds_t.lon.values
        lats = ds_t.lat.values

        hist, _, _ = np.histogram2d(lats, lons, bins=[lats_edges, lons_edges])
        ds['stats'][i] = hist

    return ds 

def mask_above_percentile(da, percentile):
    mask = da.where(da > da.quantile(percentile/100), drop=True).notnull().sum()
    # mask = da.where(da > np.percentile(da, percentile), drop=True)
    return mask

def preprocess(ds):

    ds_subset = ds.sel(latitude=slice(wesn[2], wesn[3]), longitude=slice(wesn[0], wesn[1]))
    
    return ds_subset

def get_density(ds, n_cells_deg): #, land_mask):
        
        R = 6371.0  # km
        dlat = 1/(2*n_cells_deg) # degrees(1/2 of resolution)
        dlon = 1/(2*n_cells_deg) # degrees(1/2 of resolution)

        datasets_t = []
        for t in tqdm(ds.time):
            # Define patch domain 
            ds_t = ds.sel(time=t)

            # stats is zero where land_mask is True
            # cond = land_mask['mask'].values
            # ds_t['stats'] = ds_t['stats'].where(cond, ds_t['stats'].values, 0)
            ds_t['stats'] = ds_t['stats'].where(ds_t['stats'].values, 0)

            # sum of particles in the domain
            particles = np.sum(ds_t['stats'].values)

            # Calculate area of each grid cell
            lon_grid, lat_grid = np.meshgrid(ds_t['stats'].lon, ds_t['stats'].lat) 
            area = R**2 * (np.sin(np.radians(lat_grid + dlat)) - np.sin(np.radians(lat_grid - dlat))) * (np.radians(lon_grid + dlon) - np.radians(lon_grid - dlon))  
            ds_t['area'] = (('lat', 'lon'), area)
            ds_t['mass'] =ds_t['stats'] * (7e7 / particles) # in kg
            ds_t['density'] = ds_t['mass'] /ds_t['area'] # in kg per km^2

            datasets_t.append(ds_t)

        filtered_ds = xr.concat(datasets_t, dim='time')

        return filtered_ds

def get_patch_density(ds, n_cells_deg, p): # , land_mask):
        '''
        This function calculates statistics for a patch of particles in the domain following the dichotomy method.
        Main application is the NPGP, defined as the 75th percentile of the total particles in the domain, where the domain is commented bewow. 
        '''
        
        R = 6371.0  # km
        dlat = 1/(2*n_cells_deg) # degrees(1/2 of resolution)
        dlon = 1/(2*n_cells_deg) # degrees(1/2 of resolution)
        # p = 75 # percentile
        
        #wesn = [-155, -125, 20, 45]
        # ds = ds.sel(latitude=slice(wesn[2], wesn[3]), longitude=slice(wesn[0], wesn[1]))

        datasets_t = []
        for t in tqdm(ds.time):
            # Define patch domain 
            ds_t = ds.sel(time=t)

            # stats is zero where land_mask is True
            # cond = land_mask['mask'].values
            # ds_t['stats'] = ds_t['stats'].where(cond, ds_t['stats'].values, 0)
            ds_t['stats'] = ds_t['stats'].where(ds_t['stats'].values, 0)

            particles = np.sum(ds_t['stats'].values) # Total particles in the domain
            # Return a copy of the array collapsed into one dimension
            sorted_count = np.sort(ds_t['stats'].values.flatten())[::-1] 
            cumsum = sorted_count.cumsum() 
            patch_particles = (0.75 * particles).astype(int)
            target = p/100 * particles
            idx = np.argmax(cumsum >= target)
            # Define the range of values in the patch
            ps = sorted_count[:idx] 
            
            # Create a mask to define the patch as per dichotomy method
            mask = (ds_t['stats'] >= ps.min()) 
            filtered_ds_t = ds_t.where(mask) # keep same domain, send filtered to nan
            
            # Calculate area of each grid cell
            lon_grid, lat_grid = np.meshgrid(filtered_ds_t['stats'].lon, filtered_ds_t['stats'].lat) 
            area = R**2 * (np.sin(np.radians(lat_grid + dlat)) - np.sin(np.radians(lat_grid - dlat))) * (np.radians(lon_grid + dlon) - np.radians(lon_grid - dlon))  
            filtered_ds_t['area'] = (('lat', 'lon'), area)
            filtered_ds_t['area'] = filtered_ds_t['area'].where(filtered_ds_t['stats'].notnull(), np.nan) # area to nan if stats is nan
            filtered_ds_t['mass'] = filtered_ds_t['stats'] * (7e7 / patch_particles) # in kg
            filtered_ds_t['density'] = filtered_ds_t['mass'] / filtered_ds_t['area'] # in kg per km^2
            filtered_ds_t['stats'] = filtered_ds_t['stats'] / particles # normalised particles

            datasets_t.append(filtered_ds_t)

        filtered_ds = xr.concat(datasets_t, dim='time')

        return filtered_ds