"""
advect on HYCOM surface currents
"""
import sys

sys.path.insert(1, "./src")

from src.run_advector_2D import run_advector_2D
from src.preprocessors import rename_var_preprocessor
from get_previous_advections import make_sourcefile_from_advection_output
from Mattia_GLOBCURRENT_UNIFORM import readme_txt
from datetime import datetime, timedelta
from pathlib import Path

# HYCOM_U_PATH = "I:/ADVECTOR/metocean/HYCOM_ADVECTOR_FORMATTED_1993_2020/u/u_*.nc"
# HYCOM_V_PATH = "I:/ADVECTOR/metocean/HYCOM_ADVECTOR_FORMATTED_1993_2020/v/v_*.nc"
sourcefile = "Z:/ADVECTOR/sources/uniform/T0-8km_NPO_RD2021_landmask.nc"
HYCOM_PATH = "V:/_history_reanalysis_3y/hindcast/Global_Ocean_Models/hycom/*.nc"
# advector_path = "R:/PROJECTS/TRAPs/ADVECTOR/outputs/advector_v1_2d_hycom_uniform_GPGP_1x_wind_0_2020_2021/ADVECTOR_2D_output_2020.nc"
# sourcefile = "R:/PROJECTS/TRAPs/ADVECTOR/outputs/advector_v1_2d_hycom_uniform_GPGP_1x_wind_0_2020_2021/sources_2021.nc"
# make_sourcefile_from_advection_output(Path(advector_path), Path(sourcefile))
# ADVECTION_START = datetime(2020, 12, 31)
ADVECTION_START = datetime(2021, 1, 1)
ADVECTION_END = datetime(2023, 12, 31)
WINDAGE_COEFF = 0
OUTPUTFILE_PATH = f"R:\PROJECTS\TRAPs\ADVECTOR\outputs_2\\advector_v1_2d_hycom_T0-8km_NPO_RD2021_landmask_GPGP_1x_wind_{WINDAGE_COEFF}_{ADVECTION_START.year}_{ADVECTION_END.year}"
dataset_id = None
water_varname_map = {"water_u": "U", "water_v": "V"}

if __name__ == "__main__":

    readme_content = f"""

    Description:
    GLOBCURRENT {ADVECTION_START} - {ADVECTION_END} ADVECTOR run. 

    Input:
    Located in {HYCOM_PATH}

    Sources:
    {sourcefile}

    Contributors:
    - M. Romero 
        """

    readme_txt(OUTPUTFILE_PATH, ADVECTION_START, ADVECTION_END, dataset_id, sourcefile)

    out_paths = run_advector_2D(
        output_directory=OUTPUTFILE_PATH,
        sourcefile_path=sourcefile,
        u_water_path=HYCOM_PATH,
        v_water_path=HYCOM_PATH,
        water_preprocessor=rename_var_preprocessor(water_varname_map),
        windage_coeff=WINDAGE_COEFF,
        eddy_diffusivity=0.1,  # m^2 / s
        advection_start_date=ADVECTION_START,
        timestep=timedelta(hours=1),
        num_timesteps=24 * (ADVECTION_END - ADVECTION_START).days,
        save_period=24,
        memory_utilization=0.2,  # decrease if RAM overloaded.  Can be close to 1 on dedicated compute device (e.g. GPU)
        opencl_device=None,
    )
