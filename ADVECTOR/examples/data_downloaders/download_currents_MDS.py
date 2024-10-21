import glob
import subprocess
from pathlib import Path

import pandas as pd
import xarray as xr
from tqdm import tqdm
import copernicusmarine as copernicusmarine

COPERNICUS_USER = "ypham"  # "apeytavin"
COPERNICUS_PASSWORD = "P7f6DR74iK4pW3f"  # "J2#?uQtHTxC_?an"

def download_ReproGeE_daily_surface_MDS(date: pd.Timestamp, var: str, destination: Path):
    # This script downloads a GLORYS monthly file using Copernicus Marine Toolbox and requires :
    #  - an account at https://resources.marine.copernicus.eu/
    #  - the filename and folder of the output

    dataset_id="cmems_obs_mob_glo_phy-cur_my_0.25deg_P1D-m"
    minimum_longitude = -179.875
    maximum_longitude = 179.875 
    minimum_latitude = -89.875
    maximum_latitude = 89.875
    minimum_depth = 0.
    maximum_depth = 0.
    start_datetime = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00:00'
    end_datetime = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00:00'
    variables = var, 
    output_directory = destination.parent.as_posix()
    output_filename = destination.name
    username = COPERNICUS_USER
    password = COPERNICUS_PASSWORD

    copernicusmarine.subset(dataset_id=dataset_id,
    minimum_longitude=minimum_longitude, maximum_longitude=maximum_longitude,
    minimum_latitude=minimum_latitude,maximum_latitude=maximum_latitude,
    minimum_depth=minimum_depth,maximum_depth=maximum_depth,
    start_datetime=start_datetime,end_datetime=end_datetime,
    variables=variables,
    output_directory=output_directory,output_filename=output_filename,
    username=username,password=password,
    force_download=True)#, overwrite_metadata_cache= True)


def download_and_process_currents(out_dir: Path, tstart: str, tend: str):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    while True:
        none_failed = True
        for day in pd.date_range(tstart, tend, freq="D"):
            for varname in ("uo", "vo"):
                filename = f"{varname}_{day.strftime('%Y_%m_%d')}.nc"

                # download the raw data
                raw_out_path = raw_dir / filename

                if raw_out_path.is_file():
                    try:
                        ds = xr.open_dataset(raw_out_path)
                    except (OSError, ValueError, RuntimeError, KeyError):
                        # file didn't download properly, delete and try again later
                        none_failed = False
                        raw_out_path.unlink()
                        continue
                else: 
                    ds = download_ReproGeE_daily_surface_MDS(day,varname,raw_out_path)


        if none_failed:
            print("All files downloaded.")
            break