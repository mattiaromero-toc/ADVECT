"""
Author: Mattia Romero, The Ocean Cleanup
Date: 2024-10-29
Description: Create a uniform grid of sources over the North Pacific Ocean and remove land points

Details:
- Details on the Fibonacci sphere https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf and https://devforum.roblox.com/t/generating-equidistant-points-on-a-sphere/874144 
- Optional Fibonacci sphere package https://github.com/ksbg/equidistantpoints 
- Domain coordinates of North Pacific Ocean https://www.marineregions.org/gazetteer.php/gazetteer.php?p=details&id=1908 
- The land mask is created using the global_land_mask package https://github.com/toddkarin/global-land-mask
- The sources are saved in a NetCDF file format dataframe: /Volumes/storage4/ADVECTOR/sources/uniform/T1-500_RD2020.nc  
"""

from pathlib import Path as path

import numpy as np
import pandas as pd
import xarray as xr

from global_land_mask import globe

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def fibonacci_sphere(samples=1000):
    """
    Generate `samples` points uniformly distributed on a sphere using the Fibonacci lattice.
    
    Parameters:
        samples (int): Number of points to generate.

    Returns:
        (lat, lon): Tuple of arrays with latitudes and longitudes of points in degrees.

    Conventions:
        - Latitude is the angle from the z-axis with range [-90, 90] degrees.
        - Longitude is the angle from the x-axis with range [-180, 180] degrees.
    """
    indices = np.arange(0, samples, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / samples)  # Latitude-like angle
    theta = 2 * np.pi * indices / ((1 + np.sqrt(5)) / 2)  # Longitude-like angle with golden ratio

    # Convert to Cartesian coordinates
    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    lat = np.degrees(np.arcsin(z))  # Convert z to latitude in degrees
    lon = np.degrees(np.arctan2(y, x))  # Convert x, y to longitude in degrees

    return lat, lon

def haversine(lat1, lon1, lat2, lon2, radius=6371):
    """
    Calculate the great-circle distance between two points on the Earth's surface.
    Default radius of Earth = 6371 km.
    """
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return radius * c


if __name__ == "__main__":

    RootOut = path("/Volumes/storage4/ADVECTOR/sources/uniform")

    ### Example visualization of the Fibonacci sphere

    latitudes, longitudes = fibonacci_sphere(1e3)
    globe_land_mask = globe.is_land(latitudes, longitudes)
    latitudes_masked = latitudes[~globe_land_mask]
    longitudes_masked = longitudes[~globe_land_mask]

    # Plot on cartesian map and on 3D sphere
    RootFig = RootOut / "figures"
    RootFig.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    ax = axs[0]
    ax = plt.subplot(121, projection=ccrs.PlateCarree())
    ax.scatter(longitudes, latitudes, s=1, color='b')
    ax.scatter(longitudes_masked, latitudes_masked, s=1, color='r')
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines()
    ax.set_title("Cartesian map")
    ax = axs[0]
    ax = plt.subplot(122, projection=ccrs.Orthographic(central_latitude=45, central_longitude=0))
    ax.scatter(longitudes, latitudes, s=1, transform=ccrs.PlateCarree())
    ax.scatter(longitudes_masked, latitudes_masked, s=1, color='r', transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines()
    ax.set_title("3D sphere")
    plt.savefig(RootFig / "fibonacci_sphere_uniform_release.png")
    plt.close()


    ### Creating a uniform grid of sources over the North Pacific Ocean

    wesn = list([117, -78, 0, 67]) # [W, E, S, N]
    size_domain = haversine(wesn[2], wesn[0], wesn[2], wesn[1]) * haversine(wesn[2], wesn[0], wesn[3], wesn[0]) # km2
    wesn = [wesn[i] if wesn[i] >= 0 else (wesn[i] + 360 if i < 2 else wesn[i]) for i in range(4)] # Convert in 0, 360 longitude
    n_particles = 1e6
    releasedate = "1993-01-01 00:00:00"
    print(f"Release year: {releasedate[0:4]}")

    i = 0
    lons = np.array([])
    while lons.shape[0]/size_domain < 1e-2:
        lats, lons = fibonacci_sphere(int(n_particles*1e1**i))
        lons = np.where(lons < 0, lons + 360, lons) # Convert in 0, 360 longitude
        mask = (lons >= wesn[0]) & (lons <= wesn[1]) & (lats >= wesn[2]) & (lats <= wesn[3])
        lons, lats = lons[mask], lats[mask] 
        lons = np.where(lons > 180, lons - 360, lons) # Convert back in -180, 180 longitude
        i += 1
        print(f"Ratio particles/km2: {lons.shape[0]/size_domain}")

    print(f"Number of particles in domain is: {lons.shape[0]}")
        
    globe_land_mask = globe.is_land(lats, lons)

    lons_masked = lons[~globe_land_mask]
    lats_masked = lats[~globe_land_mask]

    print(f"Number of unmasked particles in domain is: {lons_masked.shape[0]}")
    res_km = int(np.sqrt(size_domain/lons.shape[0]))
    filename = f"T0-{res_km}km_NPO_RD{releasedate[0:4]}_landmask"

    release_date = np.ones(lons_masked.shape[0])
    release_date = np.repeat(pd.to_datetime(releasedate, format="%Y-%m-%d %H:%M:%S"), lons_masked.shape[0])

    sources_xr = xr.Dataset()
    sources_xr["lon"] = xr.DataArray(lons_masked, dims=("p_id",), coords={"p_id": np.arange(lons_masked.shape[0])})
    sources_xr["lat"] = xr.DataArray(lats_masked, dims=("p_id",), coords={"p_id": np.arange(lats_masked.shape[0])})
    sources_xr["release_date"] = xr.DataArray(release_date, dims=("p_id",), coords={"p_id": np.arange(lons_masked.shape[0])})
    # sources_xr['unsd'] = xr.DataArray(np.ones(lon_grid_masked.shape[0]), dims=("p_id",), coords={"p_id": np.arange(lon_grid_masked.shape[0])})
    sources_xr.attrs["title"] = "Particle sources for ADVECTOR-type Lagrangian Dispersion Model"
    sources_xr.attrs["institution"] = "The Ocean Cleanup"
    sources_xr.attrs["source"] = f"Uniform grid with {res_km} km resolution over the North Pacific Ocean"
    sources_xr.to_netcdf(path(f"/Volumes/storage4/ADVECTOR/sources/uniform/{filename}.nc"))

    # plot 
    fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    sources_xr.plot.scatter(x="lon", y="lat", ax=ax, color='r', s=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', alpha=0.5)
    plt.savefig(RootFig / f"{filename}.png")
    plt.close()

    print(f"Sources saved in {RootOut / filename}.nc")

    print("Enjoy!")

