import numpy as np

import pdb

class InDomainError(Exception):
    pass


class BoundaryConditions:
    _dims = ['x', 'y', 'z']
    _directions = {'x': ['west', 'east'], 'y': ['south', 'north'], 'z': ['bottom', 'top']}

    def __init__(self, x_west=0, x_east=0, y_south=0, y_north=0, z_bottom=0, z_top=0):
        # Initialize the dictionary of boundary conditions
        self.boundary_cond_dict = dict()
        self.boundary_cond_dict['x_west'] = self._make_boundary_method(x_west)
        self.boundary_cond_dict['x_east'] = self._make_boundary_method(x_east)
        self.boundary_cond_dict['y_south'] = self._make_boundary_method(y_south)
        self.boundary_cond_dict['y_north'] = self._make_boundary_method(y_north)
        self.boundary_cond_dict['z_bottom'] = self._make_boundary_method(z_bottom)
        self.boundary_cond_dict['z_top'] = self._make_boundary_method(z_top)

    def get_boundary_point(self, values, indices):
        boundary_fxn = self._get_boundary_function(values.shape(), indices)
        return boundary_fxn(values, indices)

    @staticmethod
    def _make_boundary_method(value):
        # All boundary methods must take as arguments the matrix of values (in the shape of the domain)
        # and the indices of the ghost point needed as a tuple. They must return the value of that ghost point.

        # If the value given is a scalar value, just return that value (as a float)
        if callable(value):
            return value
        elif isinstance(value, (int, float)):
            return lambda c, inds: float(value)
        else:
            raise NotImplementedError('No method implemented to make a boundary condition for value of type {}'.format(type(value)))

    def _get_boundary_function(self, shape, indices):
        out_of_domain_index = -1
        # Search the indices tuple for the dimension that is out of the domain. Only one dimension may be out of the
        # domain
        for i in range(len(indices)):
            if indices[i] < 0 or indices[i] > shape[i]:
                if out_of_domain_index >= 0:
                    raise IndexError('The requested index is out of the domain in more than one dimension')
                else:
                    out_of_domain_index = i
                    # Use this to construct the proper key for the boundary condition. The class attribute
                    # _directions is a dictionary with each dimension name as a key and each value a list where
                    # the first value is the negative side direction and the second the positive side direction
                    if indices[i] < 0:
                        direction_index = 0
                    else:
                        direction_index = 1

        if out_of_domain_index < 0:
            # If this is still true, then our point must actually be in the domain.
            raise InDomainError('The point {} lies within the domain'.format(indices))

        dim = self._dims[out_of_domain_index]
        direction = self._directions[dim][direction_index]
        function_key = '{}_{}'.format(dim, direction)
        return self.boundary_cond_dict[function_key]


def check_transport_inputs(dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions):
    """
    Verify that all the inputs for any of the transport operators are the correct types and shapes
    :param dt:
    :param dx:
    :param dy:
    :param dz:
    :param u_x:
    :param u_y:
    :param u_z:
    :param D_x:
    :param D_y:
    :param D_z:
    :return: all parameters, converted to proper types as necessary (i.e. scalar numbers are forced to floats)
    """

    must_be_scalar_msg = '{} must be a scalar int or float'

    if isinstance(dt, (int, float)):
        dt = float(dt)
    else:
        raise TypeError(must_be_scalar_msg.format('dt'))

    if isinstance(dx, (int, float)):
        dx = float(dx)
    elif dx is not None:
        raise TypeError(must_be_scalar_msg.format('dx'))

    if isinstance(dy, (int, float)):
        dy = float(dy)
    elif dy is not None:
        raise TypeError(must_be_scalar_msg.format('dy'))

    if isinstance(dz, (int, float)):
        dz = float(dz)
    elif dz is not None:
        raise TypeError(must_be_scalar_msg.format('dz'))

    # Now, if u_x is a scalar, then all the u's and D's must be scalars too, and domain_size must be given with the
    # correct number of dimension lengths. If u_x is a numpy array, then all other u's and D's must be numpy arrays
    # of the same shape.
    u_and_D = [u_x, u_y, u_z, D_x, D_y, D_z]
    if isinstance(u_x, (int, float)):
        is_scalar = True
        check_types = [isinstance(el, (int, float)) or el is None for el in u_and_D]
        if not all(check_types):
            raise TypeError('If u_x is a scalar, then u_y, u_z, D_x, D_y, and D_z must all be scalar or None')
        elif domain_size is None:
            raise ValueError('If the u and D values are scalars, then domain_size must be given')
        elif not isinstance(domain_size, tuple):
            raise TypeError('domain_size must be given as a tuple')

        # Ensure all are floats, unless the optional ones are None
        u_x = float(u_x)
        D_x = float(D_x)
        u_y = float(u_y) if u_y is not None else None
        D_y = float(D_y) if D_y is not None else None
        u_z = float(u_z) if u_z is not None else None
        D_z = float(D_z) if D_z is not None else None

    elif isinstance(u_x, np.ndarray):
        is_scalar = True
        check_types = [isinstance(el, np.ndarray) or el is None for el in u_and_D]
        if not all(check_types):
            raise TypeError('If u_x is a numpy array, then u_y, u_z, D_x, D_y, and D_z must all be numpy arrays or None')

        check_sizes = [el.shape == u_x.shape for el in u_and_D if el is not None]
        if not all(check_sizes):
            raise ValueError('If u_x is a numpy array, then all u and D values must have the same shape or be None')

        # Set domain_size equal to the shape of the arrays so that the calling function can use that.
        domain_size = check_sizes[0]

        # Ensure all are 64 bit floats
        u_x = u_x.astype(np.float64)
        D_x = D_x.astype(np.float64)
        u_y = u_y.astype(np.float64) if u_y is not None else None
        D_y = D_y.astype(np.float64) if D_y is not None else None
        u_z = u_z.astype(np.float64) if u_z is not None else None
        D_z = D_z.astype(np.float64) if D_z is not None else None

    else:
        raise TypeError('u_x must be a scalar int or float, or a numpy array')

    # Check the dimensionality. Both u_x and D_x must be given. If u_y is given, then D_y must also be given. If u_z is
    # given, then D_z must also be given. If u_y/D_y are NOT given, then u_z/D_z must NOT be given. Check that the right
    # number of dimensions is given, whether that is by the number of dims in the u and D numpy arrays or in the
    # domain_size tuple.
    model_dims = 1
    if u_x is None or D_x is None:
        raise ValueError('u_x and D_x must both be given')

    if (u_y is None) ^ (D_y is None):
        raise ValueError('Both or neither u_y and D_y must be given')
    elif u_y is not None:
        model_dims = 2

    if (u_z is None) ^ (D_z is None):
        raise ValueError('Both or neither u_z and D_z must be given')
    elif u_z is not None:
        model_dims = 3

    if u_z is not None and u_y is None:
        raise ValueError('Cannot give u_z/D_z and not u_y/D_y - 2D models should use the x and y dimensions')

    if model_dims >= 2 and dy is None:
        raise ValueError('Based on the u and D values given, this is a {0}D model, so dy must be given'.format(model_dims))

    if model_dims >= 3 and dz is None:
        raise ValueError('Based on the u and D values given, this is a {0}D model, so dz must be given'.format(model_dims))

    if is_scalar:
        if len(domain_size) != model_dims:
            raise ValueError('Based on the u and D values given, this is a {0}D model, so domain_size must have {0} elements'.format(model_dims))
    else:
        check_dims = [el.ndim == model_dims for el in u_and_D if el is not None]
        if not all(check_dims):
            raise ValueError('Based on the u and D values given, this is a {0}D model, so all u and D values must have {0} dimensions'.format(model_dims))

    if boundary_conditions is None:
        # If not given, initialize boundary conditions to the default (0 values on all edges).
        boundary_conditions = BoundaryConditions()
    elif not isinstance(boundary_conditions, BoundaryConditions):
        raise TypeError('boundary_conditions must be an instance of transport.BoundaryConditions')

    # Finally check values. dt, dx, dy, dz, D_x, D_y, and D_z must all be positive, dt, dx, dy, dz must be > 0
    # u_x, u_y, and u_z may be positive or negative, so we do not need to check.
    if dt <= 0:
        raise ValueError('dt must be > 0')
    elif dx <= 0:
        raise ValueError('dx must be > 0')
    elif model_dims >= 2 and dy <= 0:
        raise ValueError('dy must be > 0')
    elif model_dims >= 3 and dz <= 0:
        raise ValueError('dz must be > 0')

    if np.any(D_x < 0):
        raise ValueError('D_x must be >= 0')
    elif model_dims >= 2 and np.any(D_y < 0):
        raise ValueError('D_y must be >= 0')
    elif model_dims >= 3 and np.any(D_z < 0):
        raise ValueError('D_z must be >= 0')

    return dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions, model_dims


def add_coefficient_to_row(row, coeff, curr_lin_idx, domain_size, di=0, dj=0, dk=0):
    try:
        row[_get_relative_linear_index(curr_lin_idx, domain_size, di=di, dj=dj, dk=dk)] += coeff
    except ValueError:
        # A ValueError occurs if we try to index a point outside the domain.
        # Skip this for now and deal with it as a boundary condition
        pass

    # Do not need to return, row is modified in place


def _get_relative_linear_index(curr_lin_idx, dims, di=0, dj=0, dk=0):
    """
    Return the linear index for a point some number of grid cells away from the current one
    :param curr_lin_idx: the current grid cell's linear index, i.e. the result of
    numpy.ravel_multi_index((i,j,k), dims=c.shape), where c is an array of the physical quantity in the shape of the
    domain
    :param dims: the dimensions of the domain, as a sequence (list, tuple, or numpy array)
    :param di: the number of grid cells to move away in the x direction (positive or negative).
    :param dj: the number of grid cells to move away in the x direction (positive or negative). dims must have at least
    two elements for this to be nonzero.
    :param dk: the number of grid cells to move away in the x direction (positive or negative). dims must have at least
    three elements for this to be nonzero.
    :return: the linear index of the grid cell (di, dj, dk) away from the current one.
    """

    if not isinstance(curr_lin_idx, int):
        raise TypeError('curr_lin_idx must be a int')
    if not isinstance(di, int):
        raise TypeError('di must be a int')
    if not isinstance(dj, int):
        raise TypeError('dj must be a int')
    if not isinstance(dk, int):
        raise TypeError('dk must be a int')

    dims = np.array(dims)
    if dims.dtype != np.int:
        raise TypeError('dims must be a sequence of ints')
    elif dims.ndim != 1:
        raise TypeError('dims must be a 1D sequence')

    if dims.size < 2 and dj != 0:
        raise TypeError('dims must have 2 or more elements for dj to be nonzero')
    if dims.size < 3 and dk != 0:
        raise TypeError('dims must have 3 or more elements for dk to be nonzero')

    # Convert the linear index for the current grid cell to multi-index notation, then add [di, dj, dk] to it to find
    # the new request point. Then reconvert back to linear index.
    curr_idx = np.array(np.unravel_index(curr_lin_idx, dims=dims, order='F'))
    dijk = [di, dj, dk]
    dijk = np.array(dijk[:dims.size])  # make dijk the same length as curr_idx
    return np.ravel_multi_index(curr_idx + dijk, dims=dims, order='F')