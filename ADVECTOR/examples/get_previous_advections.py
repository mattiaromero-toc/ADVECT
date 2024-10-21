from os import walk
import xarray as xr
from datetime import date
from pathlib import Path


def make_sourcefile_from_advection_output(
    advection_file_path: Path, sourcefile_path: Path
):
    ds_output = xr.open_dataset(advection_file_path.as_posix())
    ds_sources = ds_output.isel(time=-1).drop("time").squeeze()
    ds_sources.to_netcdf(sourcefile_path.as_posix())


def get_output_file_for_sources(
    outputs_folder: Path, dispersion_2020_file_path: Path, metocean_model: str
):
    _, _, output_files = next(walk(outputs_folder.as_posix()))
    computed_dates = sorted(
        [
            f.split(f"operational_advector_output_")[1].strip(".nc")
            for f in output_files
            if f.find(f"operational_advector_output_") >= 0
        ]
    )

    if len(computed_dates) >= 1:
        return (
            outputs_folder
            / Path("operational_advector_output_" + computed_dates[-1] + ".nc"),
            date.fromisoformat(computed_dates[-1]),
        )
    else:
        return (
            dispersion_2020_file_path,
            date(year=2020, month=12, day=31)
            if metocean_model == "hycom"
            else date(year=2019, month=12, day=31),
        )