"""
The ``chem_solver`` module contains the functions necessary to initialize and solve idealized chemistry cases.
"""
import os
import netCDF4 as ncdf
from ..utilities.config import ConfigurationError, get_domain_size_from_config
from ..chemderiv import rhs
import numpy as np
from scipy.integrate import odeint
import time
_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
mech_dir = os.path.join(_mydir, '..', 'Mechanisms')
input_dir = os.path.join(_mydir, '..', '..', 'inputs')

def _parse_pecan_species(spec_file):
    """
    Parse a PECANS-style species file. The PECANS format requires that each specie is defined on
    its own line. Anything following a # is considered a comment and ignored
    :param spec_file: The path to the species file, as a str
    :return: nothing. All species are added to Species.instances.
    """
    species_name = []
    with open(spec_file, 'r') as f:
        for line in f:
            try:
                i = line.index('#')
            except ValueError:
                line2 = line.strip()
            else:
                line2 = line[:i].strip()

            if len(line2) > 0:
                species_name.append(line2)
    return tuple(species_name)

def init_explicit_nox_chem_solver(config):
    """
    Initialization function for the simplified nox-o3 chemistry solver.

    :param kwargs: extra keyword arguments not used in this function.

    :return: the solver function and a tuple of species names.
    :rtype: function and tuple of str
    """
    species_file = os.path.join(mech_dir, config.get('CHEMISTRY', 'mechanism') + '.spc')
    species = _parse_pecan_species(species_file)
    if 'fixed_params' in config.section_as_dict('CHEMISTRY'):
        temp = config.get('CHEMISTRY', 'fixed_params')['temp']
        cair = config.get('CHEMISTRY', 'fixed_params')['nair']
    else:
        msg = 'Temperature and air densities are not provided.'
        raise ConfigurationError(msg)

    const_param = [temp, cair]
    const_species = list()
    chem_const_species = dict()
    if 'const_species' in config.section_as_dict('CHEMISTRY'):
        const_species.extend(list(config.get('CHEMISTRY', 'const_species')))
        const_param.extend(list(config.get('CHEMISTRY', 'const_species').values()))
        for specie in const_species:
            chem_const_species[specie] = np.zeros(get_domain_size_from_config(config)) + \
                                         config.get('CHEMISTRY', 'const_species')[specie]
    forced_param = []
    if 'forced_species' in config.section_as_dict('CHEMISTRY'):
        forced_species = list(config.get('CHEMISTRY', 'forced_species'))
        const_species.extend(forced_species)
        if 'forced_input' in config.section_as_dict('CHEMISTRY'):
            f = ncdf.Dataset(os.path.join(input_dir, config.get('CHEMISTRY', 'forced_input')))
            for specie in forced_species:
                forced_param.append(np.array(f.variables[specie]))
                chem_const_species[specie] = np.array(f.variables[specie])
        else:
            msg = 'The input files for forced species are not provided.'
            raise ConfigurationError(msg)
    species = tuple([specie for specie in species if specie not in const_species])
    def chem_solver(dt, const_param, forced_param, **species_in):
        dt = float(dt)
        lens = [len(conc) for species, conc in species_in.items()][0]
        species = list(species_in.keys())
        init_time = time.time()
        for i in range(lens):
            this_param = const_param + [param[i] for param in forced_param]
            this_conc = [conc[i] for species, conc in species_in.items()]
            sol = odeint(rhs, np.array(this_conc), np.arange(0, dt, 1), args=(np.array(this_param),))
            for j in range(len(species)):
                # if this_species not in const_species.keys():
                species_in[species[j]][i] = sol[-1,j]
        #print("The chemistry solver took: %.6s seconds in this time step" % (time.time() - init_time))
        return species_in

    return chem_solver, species, const_param, forced_param, chem_const_species
