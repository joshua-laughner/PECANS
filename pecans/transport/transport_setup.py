import numpy as np
import pdb

from . import transport_utils as tutils, backwards_euler, crank_nicholson
from ..utilities.config import get_domain_size_from_config


def setup_transport(config):
    transport_method = get_solver(config)

    # Also set up the function to get winds for the current time since model start. If winds are fixed, this means it
    # will just always return the matrices of fixed winds
    winds_method = _get_wind_load_method(config)

    def transport_solver(values, dt, dx=None, dy=None, dz=None, u_x=None, u_y=None, u_z=None, D_x=None, D_y=None, D_z=None, domain_size=None, boundary_conditions=None):
        # Probably what I should ultimately do is create a class that represents a domain that stores a dictionary of
        # species, as well as whether each one is subject to transport and/or chemistry alongside the domain configuration
        # It can have a method to return the configuration as a dict, or store it internally as one, that can be passed in
        # here as a **kwargs. The solver or a driver function could then loop over the species to transport.
        transport_matrix = transport_method(dt=dt, dx=dx, dy=dy, dz=dz, u_x=u_x, u_y=u_y, u_z=u_z, D_x=D_x, D_y=D_y, D_z=D_z,
                              domain_size=domain_size, boundary_conditions=boundary_conditions)
        #pdb.set_trace()
        current_concentrations = tutils.domain_to_vector(values)
        next_concentrations = np.linalg.solve(transport_matrix, current_concentrations)
        return tutils.vector_to_domain(next_concentrations, domain_size)

    return transport_solver, winds_method


def get_solver(config):
    method = config.get('TRANSPORT', 'scheme')
    if method == 'backwards_euler_2' or method == 'implicit2':
        return backwards_euler.construct_transport_matrix_with_stencil
    elif method == 'crank_nicholson':
        return crank_nicholson.construct_transport_matrix_with_stencil
    else:
        raise ValueError('No solver defined for method == "{}"'.format(method))


def _get_wind_load_method(config):
    method = config.get('TRANSPORT', 'wind_type')
    domain_size = get_domain_size_from_config(config)
    if method == 'fixed':
        winds_xyz = config.get('TRANSPORT', 'wind_speeds')
        diffusion_xyz = config.get('TRANSPORT', 'diffusion_coeffs')

        u_x = winds_xyz['x']
        D_x = diffusion_xyz['x']

        if len(domain_size) >= 2:
            u_y = winds_xyz['y']
            D_y = diffusion_xyz['y']
        else:
            u_y = None
            D_y = None

        if len(domain_size) >= 3:
            u_z = winds_xyz['z']
            D_z = diffusion_xyz['z']
        else:
            u_z = None
            D_z = None

        def constant_winds(time_since_model_start):
            return u_x, u_y, u_z, D_x, D_y, D_z

        return constant_winds
    else:
        raise NotImplementedError('No winds load method defined for option "{}"'.format(method))
