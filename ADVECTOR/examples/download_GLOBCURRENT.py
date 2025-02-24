from data_downloaders.download_currents_MDS import download_and_process_currents
from datetime import datetime, timedelta
from pathlib import Path

ADVECTION_START = datetime(2023, 1, 1)
ADVECTION_END = datetime(2025, 1, 1)

if __name__ == "__main__":

 
    import os
    folder_path = "V:/metocean/GLOBCURRENT/raw"
    prefix = ["uo_2023", "vo_2023"]

    # List all files in the folder
    for filename in os.listdir(folder_path):
        if (filename.startswith(prefix[0])) or (filename.startswith(prefix[1])):
            file_path = os.path.join(folder_path, filename)
            os.remove(file_path)  
            print(f"Deleted: {file_path}")

    print("Cleanup complete.")

    # download surface currents 
    download_path = Path("V:/metocean/GLOBCURRENT") # storageADVECTOR 
    download_and_process_currents(download_path, ADVECTION_START, ADVECTION_END - timedelta(days=1))