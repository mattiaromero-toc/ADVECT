import os
from pathlib import Path 
import json 
import xarray as xr
import numpy as np
import logging
from tqdm import tqdm
import copernicusmarine as copernicusmarine

def load_config():
    """Load configuration from a JSON file."""
    config_file = Path(input("Path to config file: "))
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            return json.load(f)
    else:
        print("The entered path does not exist.")
        quit()

class CopernicusDataDownloader:
    def __init__(self):
        self.config = load_config()
    
    def confirm_download(self):
        print("\nDownload Configuration:")
        for key, value in self.config.items():
            print(f"{key}: {value}")
        
        confirm = input("Do you want to proceed with the download? (yes/no): ").strip().lower()
        return confirm.lower() == "yes"

    def setup_logging(self):
        """Sets up a warning log file."""
        logging.basicConfig(
            filename=os.path.join(self.output_dir, "warnings.log"),
            level=logging.WARNING,
            format="%(asctime)s - %(levelname)s - %(message)s"
        )

    def download_CMT(self) -> None:
        """
        Download files using Copernicus Marine Toolbox and parameters from the config file.
        This script downloads a  and requires :
        - an account at https://resources.marine.copernicus.eu/
        """

        dataset_id = self.config["dataset_id"]
        output_dir = self.config["output_directory"]
        os.makedirs(output_dir, exist_ok=True)
        
        start_datetime = self.config["start_datetime"]
        end_datetime = self.config["end_datetime"]

        # dataset_info = copernicusmarine.describe(dataset_id)
        # data_end_datetime = dataset_info["time_coverage_end"]
        # if end_datetime > data_end_datetime:
        #     end_datetime = data_end_datetime
        #     print("Final date dataset:", end_datetime)

        if not os.path.exists(f"{output_dir}/{dataset_id}.nc") and self.confirm_download():
            copernicusmarine.subset(
                username=self.config["username"], 
                password=self.config["password"],
                dataset_id=dataset_id,
                minimum_longitude=self.config["longitude"]["min"],
                maximum_longitude=self.config["longitude"]["max"],
                minimum_latitude=self.config["latitude"]["min"],
                maximum_latitude=self.config["latitude"]["max"],
                minimum_depth=self.config["depth"]["min"],
                maximum_depth=self.config["depth"]["max"],
                start_datetime=start_datetime,
                end_datetime=end_datetime,
                variables=self.config["variables"],
                output_directory=output_dir,
                output_filename=dataset_id,
                force_download=True
            )
        elif os.path.exists(f"{output_dir}/{dataset_id}.nc"):
            print("Already downloaded!")
        else: 
            print("Download cancelled.")
        
    def split_into_daily_files(self):
        """Splits downloaded dataset into daily files per variable."""
        dataset_path = Path(f'{self.config["output_directory"]}/{self.config["dataset_id"]}.nc')
        ds = xr.open_dataset(dataset_path)
        for var in ds.data_vars:
            ds_var = ds[var]
            for date, daily_data in tqdm(ds_var.groupby("time.date")):  
                file_path = f"{self.config["output_directory"]}/{var}_{date.strftime('%Y_%m_%d')}.nc"
                daily_data.to_netcdf(file_path)

            try:
                xr.open_dataset(file_path)
            except Exception as e:
                logging.warning(f"Corrupt file detected: {file_path} - {str(e)}")
                self.broken_files.append(file_path)
                os.remove(file_path)
                print(f"{file_path} is corrupted and will be reprocessed later.")

        os.remove(dataset_path)
        if self.broken_files:
            print("\nSome files were corrupted. Check warnings.log for details.")

    def run(self):
        self.download_CMT()
        self.split_into_daily_files()
        print("Processing complete!")

if __name__ == "__main__":
    
    downloader = CopernicusDataDownloader()
    downloader.run()




