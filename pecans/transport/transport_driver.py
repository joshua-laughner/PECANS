import numpy as np
import pdb

from . import transport_utils as tutils, backwards_euler


def solve(values, method, dt, dx=None, dy=None, dz=None, u_x=None, u_y=None, u_z=None, D_x=None, D_y=None, D_z=None, domain_size=None, boundary_conditions=None):
    # Probably what I should ultimately do is create a class that represents a domain that stores a dictionary of
    # species, as well as whether each one is subject to transport and/or chemistry alongside the domain configuration
    # It can have a method to return the configuration as a dict, or store it internally as one, that can be passed in
    # here as a **kwargs. The solver or a driver function could then loop over the species to transport.
    method = get_solver(method)

    transport_matrix = method(dt=dt, dx=dx, dy=dy, dz=dz, u_x=u_x, u_y=u_y, u_z=u_z, D_x=D_x, D_y=D_y, D_z=D_z,
                          domain_size=domain_size, boundary_conditions=boundary_conditions)
    current_concentrations = tutils.domain_to_vector(values)
    next_concentrations = np.linalg.solve(transport_matrix, current_concentrations)
    return tutils.vector_to_domain(next_concentrations, domain_size)


def get_solver(method):
    if isinstance(method, str):
        if method == 'backwards_euler_2' or method == 'implicit2':
            return backwards_euler.construct_transport_matrix_with_stencil
        else:
            raise ValueError('No solver defined for method == "{}"'.format(method))
    else:
        #TODO: ensure method is callable in this case
        return method