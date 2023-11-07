from argparse import ArgumentParser
import matplotlib.pyplot as plt
import netCDF4 as ncdf
import numpy as np
from pathlib import Path
import re

from typing import Sequence


def driver(output_dir, plot_file):
    output_files = sorted(Path(output_dir).glob('pecans_output_*.nc'))
    times, variables = _load_variables(output_files)
    nvar = len(variables)
    _, axs = plt.subplots(nvar, 1, figsize=(12, 4*nvar), sharex=True)

    for ax, varname in zip(axs, variables):
        ax.plot(times, variables[varname])
        ax.set_ylabel(f'{varname} (molec/cm^3)')

    ax.set_xlabel('Seconds since model start')
    plt.savefig(plot_file, bbox_inches='tight')


def _determine_plot_variables(output_file: Path) -> Sequence[str]:
    variables = []
    with ncdf.Dataset(output_file) as ds:
        for varname in ds.variables.keys():
            if varname in {'x', 'y', 'z'}:
                continue
            elif varname.startswith('E_'):
                continue
            else:
                variables.append(varname)

    return variables



def _time_from_pecans_output_name(filename: str) -> int:
    filename = Path(filename).stem
    matches = re.search(r'pecans_output_(?P<days>\d{3})d(?P<hours>\d{2})h(?P<minutes>\d{2})m(?P<seconds>\d{2})s', filename)
    days = int(matches.group('days'))
    hours = int(matches.group('hours'))
    minutes = int(matches.group('minutes'))
    seconds = int(matches.group('seconds'))

    return seconds + (minutes * 60) + (hours * 60 * 60) + (days * 60 * 60 * 24)


def _load_variables(output_files: Sequence[Path]):
    varnames = _determine_plot_variables(output_files[0])
    times = np.asarray([_time_from_pecans_output_name(f.stem) for f in output_files])
    variables = {v: np.full_like(times, np.nan) for v in varnames}
    for ifile, file in enumerate(output_files):
        with ncdf.Dataset(file) as ds:
            for v in varnames:
                variables[v][ifile] = ds[v][:].filled(np.nan)

    return times, variables

def parse_args():
    p = ArgumentParser(description='Plot PECANS output')
    p.add_argument('output_dir', help='Path to the directory containing PECANS output files')
    p.add_argument('plot_file', help='Where to save the resulting plots. The type of output file is determined from the extension (e.g. PDF, JPG, etc.)')

    return vars(p.parse_args())


def main():
    clargs = parse_args()
    driver(**clargs)


if __name__ == '__main__':
    main()
