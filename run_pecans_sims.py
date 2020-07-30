#!/opt/anaconda3/envs/PECANS-env/bin/python
import argparse
import os

import numpy as np

from pecans.ensembles.api import EnsembleRunner
from pecans.utilities.config import load_config_file

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
config_file = os.path.join(_mydir, 'pecans_config.cfg')


def name_first_order_output_files(index, **config_opts):
    lifetime_hours = config_opts['CHEMISTRY/mechanism_opts/lifetime_seconds'] / 3600
    emissions_width_km = config_opts['EMISSIONS/emission_opts/width_x'] / 1000
    return 'pecans_ens_tau-{}h_emwidth-{}km'.format(lifetime_hours, emissions_width_km)

def name_first_order_winds_output_files(index, **config_opts):
    winds = config_opts['TRANSPORT/wind_speeds/x']
    return 'pecans_ens_windspeed_{}m_s'.format(winds)

def name_two_phases_first_order_output_files(index, **config_opts):
    first_lifetime_hours = config_opts['CHEMISTRY/mechanism_opts/first_lifetime_seconds'] / 3600
    second_lifetime_horus = config_opts['CHEMISTRY/mechanism_opts/second_lifetime_seconds'] / 3600
    first_phase_width = config_opts['CHEMISTRY/mechanism_opts/first_phase_width'] / 1000
    emissions_width_km = config_opts['EMISSIONS/emission_opts/width_x'] / 1000
    return 'pecans_ens_first_tau-{}h_second_tau-{}h_fpwidth-{}km_emwidth-{}km'.format(first_lifetime_hours,
                                                                                      second_lifetime_horus,
                                                                                      first_phase_width,
                                                                                      emissions_width_km)


def sims_first_order_run_winds():
    # We want lifetimes that vary from 1-9 hours. This covers about the most extreme values we'd expect for summer NOx
    # lifetime
    winds= np.arange(3, 11, 1)

    ens = EnsembleRunner(config_file,
                         ensemble_variables={'TRANSPORT/wind_speeds/x': winds},
                         ensemble_mode='combinations',
                         save_in_individual_dirs=False,
                         save_final_output_only=True,
                         member_naming_fxn=name_first_order_winds_output_files,
                         root_output_dir=os.path.join(_mydir, '../../MATLAB/PAN_Data', 'Workspaces', 'PECANS',
                                                      'lifetime-ensemble'))

    ens.run()

def sims_first_order_run():

    # We want lifetimes that vary from 1-9 hours. This covers about the most extreme values we'd expect for summer NOx
    # lifetime
    taus = np.arange(3600, 9*3600+1, 3600)

    # We also want to test what happens when emissions widths are similar or greater than lifetimes. So we'll calculate
    # emissions widths equal to each expected lifetime
    config = load_config_file(config_file)
    winds = config.get('TRANSPORT', 'wind_speeds')
    x_wind = winds['x']
    widths = taus * x_wind
    widths = np.concatenate(([3000], widths))  # add a smaller width as an extra test

    ens = EnsembleRunner(config_file,
                         ensemble_variables={'CHEMISTRY/mechanism_opts/lifetime_seconds': taus,
                                             'EMISSIONS/emission_opts/width_x': widths},
                         ensemble_mode='combinations',
                         save_in_individual_dirs=False,
                         save_final_output_only=True,
                         member_naming_fxn=name_first_order_output_files,
                         root_output_dir=os.path.join(_mydir, '../../MATLAB/PAN_Data', 'Workspaces', 'PECANS',
                                                      'lifetime-ensemble'))

    ens.run()

def sims_two_phases_first_order_run():
    first_tau = np.arange(3600, 9*3600, 3600)
    second_tau = np.arange(3600, 9*3600, 3600)
    first_phase_width = np.arange(20*1000, 100*1000, 10*1000)
    config = load_config_file(config_file)
    winds = config.get('TRANSPORT', 'wind_speeds')
    x_wind = winds['x']
    widths = first_tau * x_wind
    widths = np.concatenate(([3000], widths))  # add a smaller width as an extra test
    #widths = [20000, 30000]

    ens = EnsembleRunner(config_file,
                         ensemble_variables={'CHEMISTRY/mechanism_opts/first_lifetime_seconds': first_tau,
                                             'CHEMISTRY/mechanism_opts/second_lifetime_seconds': second_tau,
                                             'CHEMISTRY/mechanism_opts/first_phase_width': first_phase_width,
                                             'EMISSIONS/emission_opts/width_x': widths},
                         ensemble_mode='combinations',
                         save_in_individual_dirs=False,
                         save_final_output_only=True,
                         member_naming_fxn=name_two_phases_first_order_output_files,
                         root_output_dir=os.path.join(_mydir, '../../MATLAB/PAN_Data', 'Workspaces', 'PECANS',
                                                      'lifetime-ensemble-twophases'))

    ens.run()

def main():
    parser = argparse.ArgumentParser(description='Choose one of the chemical solvers')
    parser.add_argument('solver', type=str, help='What the chemical solver is. Default is "first_order"')
    args = parser.parse_args()
    if args.solver == 'first_order':
        sims_first_order_run_winds()
    elif args.solver == 'two_phases_first_order':
        sims_two_phases_first_order_run()
    else:
        print("The chemical solver is not implemented.")
        quit()

if __name__ == '__main__':
    main()