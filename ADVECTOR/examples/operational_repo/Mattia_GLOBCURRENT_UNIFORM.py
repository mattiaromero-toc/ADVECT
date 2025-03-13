"""
ADVECTOR G+E Uniform 1993-2024
- Uniform release 1993
- Daily advection on cmems_obs_mob_glo_phy-cur_my_0.25deg_P1D G+E
- No windage 

"""
import os
import sys

sys.path.insert(1, "./src")

# from data_downloaders.download_currents_MDS import download_and_process_currents
from run_advector_2D import run_advector_2D
from preprocessors import rename_var_preprocessor
from get_previous_advections import (
    get_output_file_for_sources,
    make_sourcefile_from_advection_output,
)
from datetime import datetime, timedelta
from pathlib import Path


def readme_txt(
    out_dir,
    ADVECTION_START,
    ADVECTION_END,
    dataset_id,
    source_file_path,
    readme_content,
) -> None:

    os.makedirs(OUTPUTFILE_PATH, exist_ok=True)

    # Save the README file in the output directory
    readme_path = os.path.join(out_dir, "README.txt")
    with open(readme_path, "w") as f:
        f.write(readme_content)

    # Verify if file was created
    if os.path.exists(readme_path):
        print(f"README.txt successfully created at: {readme_path}")
    else:
        print("Error: README.txt was not created!")


ADVECTION_START = datetime(2024, 5, 31)
ADVECTION_END = datetime(2025, 1, 1)
WINDAGE_COEFF = 0
source_file_path = (
    f"I:\\ADVECTOR\\sources\\advector_studies\\T0-8km_NPO_RD1993_landmask_T1993.nc"
)
# OUTPUTFILE_PATH = f"R:\PROJECTS\\ADVECTOR Studies\\ADVECTOR\\outputs\\advector_v1_2d_globcurrent_{str(Path(source_file_path).name).split('.nc')[0]}_1x_wind_{WINDAGE_COEFF}_{ADVECTION_START.year}_{ADVECTION_END.year}"
dataset_id = [
    "cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m",
    "cmems_obs-mob_glo_phy-cur_nrt_0.25deg_P1D-m",
]
water_varname_map = {"uo": "U", "vo": "V", "longitude": "lon", "latitude": "lat"}

# K = "scale_k0.5"
# f"I:\\ADVECTOR\\sources\\advector_studies\\T0-8km_NPO_RD1993_landmask_T2003_{K}.nc"
# OUTPUTFILE_PATH = f"R:\\PROJECTS\\ADVECTOR Studies\\ADVECTOR\\outputs\\advector_v1_2d_globcurrent_T0-8km_NPO_RD1993_landmask_1x_wind_0_1992_2024\ICs\\advector_v1_2d_globcurrent_T0-8km_NPO_RD1993_landmask_T2003_{K}_1x_wind_{WINDAGE_COEFF}_{ADVECTION_START.year}_{ADVECTION_END.year}"

source_file_path = (
    f"I:\\ADVECTOR\\sources\\advector_studies\\T0-8km_NPO_RD1993_landmask_T1993.nc"
)
OUTPUTFILE_PATH = f"R:\PROJECTS\ADVECTOR Studies\ADVECTOR\outputs\\advector_v1_2d_globcurrent_T0-8km_NPO_RD1993_landmask_1x_wind_0_1992_2024"

if __name__ == "__main__":

    # pth = Path(OUTPUTFILE_PATH)
    # if not pth.exists():
    # pth.mkdir()

    # download surface currents
    download_path = f"T:\\metocean\\GLOBCURRENT\\"
    # download_and_process_currents(download_path, ADVECTION_START + timedelta(days=1), ADVECTION_END - timedelta(days=1))
    GLOBCURRENT_U_PATH = (
        download_path + "raw\\uo_*"
    )  # r"Y:\PROJECTS\TRAPs\ADVECTOR\inputs\raw\uo_*"
    GLOBCURRENT_V_PATH = (
        download_path + "raw\\vo_*"
    )  # r"Y:\PROJECTS\TRAPs\ADVECTOR\inputs\raw\vo_*"

    print(GLOBCURRENT_U_PATH)
    readme_content = f"""

    Description:
    GLOBCURRENT {ADVECTION_START} - {ADVECTION_END} ADVECTOR run. 

    Input:
    - 1993-2023: {dataset_id[0]} located in {download_path}
    - 2023-2024: {dataset_id[1]} located in {download_path}

    Sources:
    {source_file_path}

    Contributors:
    - M. Romero 
        """

    readme_txt(
        OUTPUTFILE_PATH,
        ADVECTION_START,
        ADVECTION_END,
        dataset_id,
        source_file_path,
        readme_content,
    )

    file_name = "ADVECTOR_2D_output_2024.nc"
    source_file_path = Path("I:\\ADVECTOR\\sources\\advector_studies\\") / file_name
    make_sourcefile_from_advection_output(
        Path(OUTPUTFILE_PATH) / file_name,
        source_file_path,
    )

    out_paths = run_advector_2D(
        output_directory=OUTPUTFILE_PATH,
        sourcefile_path=str(source_file_path),
        u_water_path=GLOBCURRENT_U_PATH,
        v_water_path=GLOBCURRENT_V_PATH,
        water_preprocessor=rename_var_preprocessor(water_varname_map),
        windage_coeff=WINDAGE_COEFF,
        eddy_diffusivity=0.1,  # m^2 / s
        advection_start_date=ADVECTION_START,
        timestep=timedelta(hours=1),
        num_timesteps=24 * (ADVECTION_END - ADVECTION_START).days,
        save_period=24,
        memory_utilization=0.2,  # decrease if RAM overloaded.  Can be close to 1 on dedicated compute device (e.g. GPU)
        opencl_device=None,
        show_progress_bar=True,
    )
