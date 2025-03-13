import os 
from enum import Enum
from urllib.request import urlretrieve
from urllib.error import URLError
import pandas as pd
from os import path
from pathlib import Path
from socket import timeout
from datetime import datetime
from download_currents_CMT import load_config, DataDownloader

# FIXME: https://github.com/gjpelletier/get_hycom/tree/main

class HYCOMExperiment(Enum):
    """Enumeration for different HYCOM experiments."""
    GLBy = "GLBy"
    ESPCDV02 = "ESPC-D-V02"

class HYCOMConfig():
    """Class to store and manage HYCOM experiment configurations."""

    EXPERIMENTS_STR = {
        HYCOMExperiment.GLBy: "GLBy",
        HYCOMExperiment.ESPCDV02: "ESPC-D-V02"
    }

    SUB_EXPERIMENT_RANGES = {
        HYCOMExperiment.GLBy: {
            "expt_56.3": pd.date_range(start="2014-07-01", end="2016-04-30", freq="3h"),
            "expt_57.2": pd.date_range(start="2016-05-01", end="2017-01-31", freq="3h"),
            "expt_92.8": pd.date_range(start="2017-02-01", end="2017-05-31", freq="3h"),
            "expt_57.7": pd.date_range(start="2017-06-01", end="2017-09-30", freq="3h"),
            "expt_92.9": pd.date_range(start="2017-10-01", end="2017-12-31", freq="3h"),
            "expt_93.0": pd.date_range(start="2018-01-01", end="2020-02-18", freq="3h"),
        },

        HYCOMExperiment.ESPCDV02: {
            "": pd.date_range(
                start="2024-08-10T00", end=pd.Timestamp.today().strftime("%Y-%m-%dT%H"), freq="3h"
            ),
        },
    }

    EXPERIMENTS_CONVENTION = {
        HYCOMExperiment.GLBy: (0, 360), 
        # HYCOMExperiment.GLBu: {"(-180, 180)"},
        # HYCOMExperiment.GLBv: {"(0, 360)"},
        HYCOMExperiment.ESPCDV02: (0, 360)
    }

    # Date threshold for ESPC-D-V02 experiment shift
    HYCOM_SHIFT = next(iter(SUB_EXPERIMENT_RANGES[HYCOMExperiment.ESPCDV02].values()))[0] # .strftime("%Y-%m-%dT%H")

    @classmethod
    def get_experiment_name(cls, experiment: HYCOMExperiment) -> str:
        """Returns the string name of an experiment."""
        return cls.EXPERIMENTS_STR.get(experiment, "Unknown Experiment")

    @classmethod
    def get_sub_experiment_by_time(cls, experiment: HYCOMExperiment, time: pd.Timestamp) -> str:
        """Finds the appropriate sub-experiment for a given time."""
        if experiment in cls.SUB_EXPERIMENT_RANGES:
            for sub_expt, time_range in cls.SUB_EXPERIMENT_RANGES[experiment].items():
                if time in time_range:
                    return sub_expt
        return "Unknown"

    @classmethod
    def get_experiment_convention(cls, experiment: HYCOMExperiment) -> str:
        """Returns the convention of an experiment."""
        return cls.EXPERIMENTS_CONVENTION.get(experiment, "Unknown Experiment")

class HYCOMdownload: 
    def __init__(self):
        # Initialize DataDownloader and run its process
        self.data_downloader = DataDownloader()
        self.data_downloader.run()  # This will start the download process and logging

    @classmethod
    def adjust_longitude(cls, lon, lon_range):
        """Ensures longitude values match the required convention."""
        if lon_range == (0, 360) and lon < 0:
            return lon + 360  # Convert -180:180 to 0:360
        if lon_range == (-180, 180) and lon > 180:
            return lon - 360  # Convert 0:360 to -180:180
        return lon

    @classmethod
    def make_hycom_url(cls, experiment: HYCOMExperiment, var: str, time: pd.Timestamp, wesn: list):
        """
        Generates a HYCOM daily download URL with automatic lat/lon convention adjustment.

        Args:
            experiment (HYCOMExperiment): The HYCOM experiment to use.
            time (pd.Timestamp): The timestamp of the requested data.
            three_d (bool): Whether to request 3D data.
            min_lat (float): Minimum latitude.
            max_lat (float): Maximum latitude.
            min_lon (float): Minimum longitude.
            max_lon (float): Maximum longitude.

        Returns:
            str | list[str]: The download URL(s).
        """

        sub_experiment = HYCOMConfig.get_sub_experiment_by_time(experiment, time) 
        base_url = "https://ncss.hycom.org/thredds/ncss"

        convention = HYCOMConfig.get_experiment_convention(experiment)
        if not convention:
            raise ValueError(f"Unknown experiment: {experiment}")
        
        # Adjust longitude range based on experiment convention
        wesn[0] = HYCOMdownload.adjust_longitude(wesn[0], convention)
        wesn[1] = HYCOMdownload.adjust_longitude(wesn[1], convention)

        params = (
            f"uv3z/{time.year}?var={var}&north={wesn[3]}&south={wesn[2]}&west={wesn[0]}&east={wesn[1]}"
            "&disableProjSubset=on&horizStride=1"
            f"&time_start={time.strftime('%Y-%m-%dT%H')}%3A00%3A00Z&time_end={(time + pd.Timedelta('1D') - pd.Timedelta('1s')).strftime('%Y-%m-%dT%H')}&timeStride=1&vertCoord=0&accept=netcdf4"
        ) # 2D 
        
        if experiment in [HYCOMExperiment.GLBy]: # HYCOMExperiment.GLBu, HYCOMExperiment.GLBv
            return f"{base_url}/{experiment.name}0.08/{sub_experiment}/{params}"
                    
        elif experiment == HYCOMExperiment.ESPC:
            return f"{base_url}/{experiment.name}/{params}"

    def download(self) -> None:

        output_dir = self.data_downloader.config["output_directory"]
        os.makedirs(output_dir, exist_ok=True)
        vars = self.data_downloader.config["variables"]
        start = self.data_downloader.config["start_datetime"]
        end = self.data_downloader.config["end_datetime"]
        wesn = self.data_downloader.config["wesn"]

        for time in pd.date_range(start=start, end=end, inclusive="left", freq="1D"): # NOTE: original was 6H 
            time: pd.Timestamp = time
            if time < HYCOMConfig.HYCOM_SHIFT:
                exp_enum = HYCOMExperiment.GLBy # Convert string to Enum
            else: 
                exp_enum = HYCOMExperiment.ESPCDV02

            for var in vars: 
                url = HYCOMdownload.make_hycom_url(exp_enum, var, time, wesn)
                filename = f'{var}_{time.strftime("%Y-%m-%d")}.nc'
                target_file_path = Path(f"{output_dir}/{filename}")

                if not path.exists(target_file_path):
                    print(f"Downloading {filename}")
                    # wget.download(url, str(target_file_path))

                    counter = 1
                    got_file = False
                    while (counter <= 10) and (got_file == False):
                        print('  Attempting to get data, counter = ' + str(counter))
                        tt0 = datetime.now()
                        try:
                            (a,b) = urlretrieve(url, str(target_file_path))
                            # a is the output file name
                            # b is a message you can see with b.as_string()
                        except URLError as ee:
                            if hasattr(ee, 'reason'):
                                print('  *We failed to reach a server.')
                                print('  -Reason: ', ee.reason)
                            elif hasattr(ee, 'code'):
                                print('  *The server could not fulfill the request.')
                                print('  -Error code: ', ee.code)
                        except timeout:
                            print('  *Socket timed out')
                        else:
                            got_file = True
                            print('  Downloaded data')
                        print('  Time elapsed: %0.1f seconds' % (datetime.now() - tt0).total_seconds())
                        counter += 1
                else:
                    print(f"Skipping {filename}, already exists")
                
                self.data_downloader.netCDF_qc(target_file_path)

        return print("Date range completed.")
    
    def run(self):
        self.download()
        print("Processing complete!")


if __name__ == "__main__":

    downloader = HYCOMdownload()
    downloader.run()
