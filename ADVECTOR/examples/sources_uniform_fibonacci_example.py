"""
Author: Mattia Romero, The Ocean Cleanup
Date: 2024-10-29
"""

from pathlib import Path 

import numpy as np
import pandas as pd
import xarray as xr

from global_land_mask import globe

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def fibonacci_sphere(n: int):
    """
    Generate `n` points uniformly distributed on a sphere using the Fibonacci lattice.
    
    Parameters:
        n (int): Number of points to generate.

    Returns:
        (lat, lon): Tuple of arrays with latitudes and longitudes of points in degrees.

    Conventions:
        - Latitude is the angle from the z-axis with range [-90, 90] degrees.
        - Longitude is the angle from the x-axis with range [-180, 180] degrees.
    """
    indices = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n)  # Latitude-like angle
    theta = 2 * np.pi * indices / ((1 + np.sqrt(5)) / 2)  # Longitude-like angle with golden ratio

    # Convert to Cartesian coordinates
    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    lat = np.degrees(np.arcsin(z))  
    lon = np.degrees(np.arctan2(y, x)) 

    return lat, lon

def plot_fibonacci_sphere(RootOut: Path, lons: np.array, lats: np.array) -> None:
    """
    Plot points on cartesian map and on a sphere.
    """
    
    RootFig = RootOut / "figures"
    RootFig.mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax = axs[0]
    ax.scatter(lons, lats, s=1, color='r', transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines()
    ax.set_title("Cartesian Map")
    axs[1].remove()  
    ax2 = fig.add_subplot(122, projection=ccrs.Orthographic(central_latitude=45, central_longitude=0))
    ax2.scatter(lons, lats, s=1, color='r', transform=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND)
    ax2.add_feature(cfeature.OCEAN)
    ax2.gridlines()
    ax2.set_title("Orthographic Projection")

    plt.savefig(RootFig / "fibonacci_sphere_uniform_release.png")
    plt.close()

def generate_source_file(RootOut) -> None: 
    """
    Create a uniform distribution of particles over the globe (sphere) and mask the ones over land.

    Details:
        - Fibonacci sphere method https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf 
        - The land mask is created using the global_land_mask package https://github.com/toddkarin/global-land-mask 
        - The sources are saved in a NetCDF file format dataframe.

    Returns:
        sources_xr: xr.Dataset source file for ADVECTOR. It contains a set of particles with lon, lat and release_date.  
    """
    RootOut = Path(RootOut)

    # Generate particle lons, lats 
    n = 1e4
    lats, lons = fibonacci_sphere(n)
    globe_land_mask = globe.is_land(lats, lons)
    lons_masked = lons[~globe_land_mask]
    lats_masked = lats[~globe_land_mask]

    # Visualise 
    plot_fibonacci_sphere(RootOut, lons_masked, lats_masked)

    # Generate source file 
    date = "2020-01-01 00:00:00"
    release_date = np.ones(lons_masked.shape[0])
    release_date = np.repeat(pd.to_datetime(date, format="%Y-%m-%d %H:%M:%S"), lons_masked.shape[0])
    str_release_date = release_date[0].strftime("%Y%m%d")

    sources_xr = xr.Dataset()
    sources_xr["lon"] = xr.DataArray(lons_masked, dims=("p_id",), coords={"p_id": np.arange(lons_masked.shape[0])})
    sources_xr["lat"] = xr.DataArray(lats_masked, dims=("p_id",), coords={"p_id": np.arange(lats_masked.shape[0])})
    sources_xr["release_date"] = xr.DataArray(release_date, dims=("p_id",), coords={"p_id": np.arange(lons_masked.shape[0])})

    sources_xr.attrs["title"] = "Particle sources for ADVECTOR-type Lagrangian Dispersion Model"
    sources_xr.attrs["institution"] = "The Ocean Cleanup"
    sources_xr.attrs["source"] = f"Global uniform particle release on {release_date}"

    filename = f"T{str_release_date}_global_uniform_landmask"
    filepath = RootOut/f"{filename}.nc"
    sources_xr.to_netcdf(filepath)
    print(f"Source file stored successfully as {filepath}")

if __name__ == "__main__":

    RootOut = "" # FIXME: add path 
    generate_source_file(RootOut)
