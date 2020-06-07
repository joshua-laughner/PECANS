"""
The ``chem_solver`` module contains the functions necessary to initialize and solve idealized chemistry cases.
"""
import os
from ..utilities.config import ConfigurationError
from .. import chemderiv
import numpy as np

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
mech_dir = os.path.join(_mydir, '..', 'Mechanisms')

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

def init_explicit_nox_chem_solver(mechanism='nox'):
    """
    Initialization function for the simplified nox-o3 chemistry solver.

    :param kwargs: extra keyword arguments not used in this function.

    :return: the solver function and a tuple of species names.
    :rtype: function and tuple of str
    """
    species_file = os.path.join(mech_dir, mechanism + '.spc')
    species = _parse_pecan_species(species_file)
    def chem_solver(dt, TEMP, CAIR, const_species, **species_in):
        dt = float(dt)
        lens = [len(conc) for species, conc in species_in.items()][0]
        for const_specie, const_conc in const_species.items():
            species_in[const_specie] = np.zeros(lens,) + const_conc
        for i in range(lens):
            this_conc = [conc[i] for species, conc in species_in.items()]
            this_species_in = chemderiv.chem_solver_ode(dt, TEMP, CAIR, *this_conc)
            for this_species, this_conc in this_species_in.items():
                if this_species not in const_species.keys():
                    species_in[this_species][i] = this_conc
        return species_in
    return chem_solver, species
