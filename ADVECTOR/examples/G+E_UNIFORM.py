"""
ADVECTOR G+E Uniform 2024-2024
- Uniform release 2020
- Daily advection on cmems_obs_mob_glo_phy-cur_my_0.25deg_P1D G+E
- No windage 

"""

from ADVECTOR.run_advector_2D import run_advector_2D
from data_downloaders.download_currents_MDS import download_and_process_currents
from preprocessors import rename_var_preprocessor
from get_previous_advections import (
    get_output_file_for_sources,
    make_sourcefile_from_advection_output,
)
from datetime import datetime, timedelta
from pathlib import Path


ADVECTION_START = datetime(1992, 12, 31)
ADVECTION_END = datetime(2024, 1, 1)
WINDAGE_COEFF = 0
source_file_path = r"I:\ADVECTOR\sources\uniform\T0-8km_NPO_RD2020_landmask.nc"
OUTPUTFILE_PATH = f"R:\PROJECTS\TRAPs\ADVECTOR\outputs"
OUTPUTFILE_PATH = f"R:\PROJECTS\TRAPs\ADVECTOR\outputs\\advector_v1_2d_globcurrent_T0-8km_NPO_RD2020_landmask_1x_wind_{WINDAGE_COEFF}_{ADVECTION_START.year}_{ADVECTION_END.year}"
water_varname_map = {"uo": "U", "vo": "V", "longitude": "lon", "latitude": "lat"}

if __name__ == "__main__":

    # pth = Path(OUTPUTFILE_PATH)
    # if not pth.exists():
       # pth.mkdir()

    # download surface currents 
    download_path = r"Y:\PROJECTS\TRAPs\ADVECTOR\inputs"
    download_and_process_currents(download_path, ADVECTION_START + timedelta(days=1), ADVECTION_END - timedelta(days=1))
    GLOBCURRENT_U_PATH = r"Y:\PROJECTS\TRAPs\ADVECTOR\inputs\raw\uo_*"
    GLOBCURRENT_V_PATH = r"Y:\PROJECTS\TRAPs\ADVECTOR\inputs\raw\vo_*"

    out_paths = run_advector_2D(
        output_directory=OUTPUTFILE_PATH,
        sourcefile_path=source_file_path,
        u_water_path=GLOBCURRENT_U_PATH,  
        v_water_path=GLOBCURRENT_V_PATH,  
        water_preprocessor=rename_var_preprocessor(water_varname_map),
        windage_coeff=WINDAGE_COEFF,
        eddy_diffusivity=0.1,  # m^2 / s
        advection_start_date=ADVECTION_START,
        timestep=timedelta(hours=1),
        num_timesteps=24 * (ADVECTION_END - ADVECTION_START).days,
        save_period=24,
        memory_utilization=0.4,  # decrease if RAM overloaded.  Can be close to 1 on dedicated compute device (e.g. GPU)
        opencl_device=None,
        show_progress_bar= True,
    )
