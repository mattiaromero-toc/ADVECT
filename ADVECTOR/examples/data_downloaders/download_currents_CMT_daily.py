import glob
import subprocess
from pathlib import Path

import pandas as pd
import xarray as xr
from tqdm import tqdm
import copernicusmarine as copernicusmarine

from multiprocessing import Process, Queue
from itertools import count

COPERNICUS_USER = "ypham"  # "apeytavin"
COPERNICUS_PASSWORD = "P7f6DR74iK4pW3f"  # "J2#?uQtHTxC_?an"

def download_GLOBCURRENT_daily_surface_MDS(
    date: pd.Timestamp,
    vars: list,
    tres: str,
    destination: Path
) -> None:
    # This script downloads a GLOBCURRENT daily file using Copernicus Marine Toolbox and requires :
    # - an account at https://resources.marine.copernicus.eu/
    # - var and tres
    # - the filename and folder of the output

    if tres == "H":
        tres = "PT1H-i"
    elif tres == "D":
        tres = "P1D-m"

    dataset_ids = [f"cmems_obs-mob_glo_phy-cur_my_0.25deg_{tres}", f"cmems_obs-mob_glo_phy-cur_nrt_0.25deg_{tres}"]
    minimum_longitude = -179.875
    maximum_longitude = 179.875 
    minimum_latitude = -89.875
    maximum_latitude = 89.875
    minimum_depth = 0.
    maximum_depth = 0.
    start_datetime = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00:00'
    end_datetime = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00:00'

    # FIXME: automatically extract end date 
    # dataset_info = copernicusmarine.describe(dataset_ids[0]) 
    data_end_datetime = pd.to_datetime("2023-12-31")

    if pd.to_datetime(end_datetime) > data_end_datetime and pd.to_datetime(start_datetime) < data_end_datetime:
        subsets = [
            (dataset_ids[0], start_datetime, data_end_datetime),
            (dataset_ids[1], data_end_datetime + pd.Timedelta(days=1), end_datetime),
        ]
    elif pd.to_datetime(end_datetime) > data_end_datetime and pd.to_datetime(start_datetime) > data_end_datetime:
        subsets = [(dataset_ids[1], start_datetime, end_datetime)]
    else:
        subsets = [(dataset_ids[0], start_datetime, end_datetime)]

    for ds_id, start, end in subsets:
        copernicusmarine.subset(
            dataset_id=ds_id,
            minimum_longitude=minimum_longitude, maximum_longitude=maximum_longitude,
            minimum_latitude=minimum_latitude, maximum_latitude=maximum_latitude,
            minimum_depth=minimum_depth, maximum_depth=maximum_depth,
            start_datetime=start, end_datetime=end,
            variables=[vars],  # Ensure it's a list
            output_directory=destination.parent.as_posix(),
            output_filename=destination.name,
            username=COPERNICUS_USER, password=COPERNICUS_PASSWORD,
            force_download=True,
        ) #, overwrite_metadata_cache= True)

def download_GLORYS_reanalysis_daily_surface_MDS(
    date: pd.Timestamp,
    destination: Path
) -> None:
    # This script downloads a GLORYS reanalysis daily file using Copernicus Marine Toolbox and requires :
    # - an account at https://resources.marine.copernicus.eu/
    # - var and dataset id 
    # - the filename and folder of the output

    dataset_ids = ["cmems_mod_glo_phy_my_0.083deg_P1D-m", "cmems_mod_glo_phy_myint_0.083deg_P1D-m"]
    minimum_longitude = -179.875
    maximum_longitude = 179.875
    minimum_latitude = -89.875
    maximum_latitude = 89.875
    minimum_depth = 0.
    maximum_depth = 0.
    start_datetime = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00:00'
    end_datetime = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00:00'

    # FIXME: automatically extract end date 
    # dataset_info = copernicusmarine.describe(dataset_ids[0]) 
    data_end_datetime = pd.to_datetime("2021-06-30")

    if pd.to_datetime(end_datetime) > data_end_datetime and pd.to_datetime(start_datetime) < data_end_datetime:
        subsets = [
            (dataset_ids[0], start_datetime, data_end_datetime),
            (dataset_ids[1], data_end_datetime + pd.Timedelta(days=1), end_datetime),
        ]
    elif pd.to_datetime(end_datetime) > data_end_datetime and pd.to_datetime(start_datetime) > data_end_datetime:
        subsets = [(dataset_ids[1], start_datetime, end_datetime)]
    else:
        subsets = [(dataset_ids[0], start_datetime, end_datetime)]

    for ds_id, start, end in subsets:
        copernicusmarine.subset(
            dataset_id=ds_id,
            minimum_longitude=minimum_longitude, maximum_longitude=maximum_longitude,
            minimum_latitude=minimum_latitude, maximum_latitude=maximum_latitude,
            minimum_depth=minimum_depth, maximum_depth=maximum_depth,
            start_datetime=start, end_datetime=end,
            variables=["uo","vo"],  # Ensure it's a list
            output_directory=destination.parent.as_posix(),
            output_filename=destination.name,
            username=COPERNICUS_USER, password=COPERNICUS_PASSWORD,
            force_download=True,
        ) #, overwrite_metadata_cache= True)

def download_data_with_timeout(day, process, args: list, raw_out_path) -> None:
    retry_count = 0
    counter = count(0)
    max_retries = 5   
    while retry_count < max_retries:
        p = Process(target=process, args=(day, args, raw_out_path), name='Data download')
        p.start()
        p.join(timeout=60)
        if p.is_alive():
            print("Timed out, relaunching...")
            p.terminate()
            p.join()  # Ensure it is fully terminated
            retry_count = next(counter)
        else:
            print("\nDownload completed successfully.\n")
            break
    if retry_count >= max_retries:
        raise RuntimeError("Failed to complete download after max retries.")  

def download_and_process_currents(out_dir: Path, process, args: list, tstart: str, tend: str) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    if len(args) > 1:
        raw_dir = out_dir / f"1{args[1]}_{args[2][0]}"
    elif (len(args) > 0) & (len(args) < 2):
        raw_dir = out_dir / f"1{args[1]}"
    else: 
        raw_dir = out_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    for day in pd.date_range(tstart, tend, freq="D"):
        for varname in (vars):
            filename = f"{varname}_{day.strftime('%Y_%m_%d')}.nc"
            raw_out_path = raw_dir / filename

            retry_count = 0
            counter = count(0)
            max_retries = 5   
            while retry_count < max_retries: 
                if not raw_out_path.is_file():
                    download_data_with_timeout(day, process, args, raw_out_path)
                try:
                    xr.open_dataset(raw_out_path)
                    print(f"\n{raw_out_path} passed the QC...\n")
                    break
                except Exception as e: # except (OSError, ValueError, RuntimeError, KeyError):
                    # file didn't download properly, delete and try again later
                    print(f"{raw_out_path} is corrupted, re-downloading...")
                    raw_out_path.unlink()
                    retry_count = next(counter)
                    continue

            if retry_count >= max_retries:
                raise RuntimeError("Failed to load successfully after max retries.")  


