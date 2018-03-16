import itertools
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


class StencilPoint:
    _dim_names = 'txyz'

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, value):
        if not isinstance(value, tuple) or len(value) < 2 or len(value) > 4 or any([not isinstance(el, int) for el in value]):
            raise TypeError('indices must be a 2 to 4 element tuple of ints')
        self._indices = self.expand_indices_to_dim(value)

    @property
    def coefficient(self):
        return self._coefficient

    @coefficient.setter
    def coefficient(self, value):
        if isinstance(value, int):
            value = float(value)
        elif not isinstance(value, float):
            raise TypeError('coefficient must be a scalar int or float')

        self._coefficient = value

    def __init__(self, indices, coefficient=1.0, dim=None):
        self._indices = (0, 0)
        self._coefficient = 0.0

        self.indices = self.expand_indices_to_dim(indices, dim=dim)
        self.coefficient = coefficient

    def __eq__(self, other):
        """
        Returns whether two instances of StencilPoint are equal or not. Two instances are equal if their indices are the
        same
        :param other: the other object to compare. If not a StencilPoint, then this method will always return false.
        :return: bool
        """
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.indices == other.indices

    def __str__(self):
        return '{} at {}'.format(self.coefficient, self.indices)

    def get_nonzero_dims(self, include_time=True):
        first_idx = 0 if include_time else 1
        nz_dims = [self._dim_names[i] for i in range(first_idx, len(self.indices)) if self.indices[i] != 0]
        return ''.join(nz_dims)

    def get_dim_index(self, dim_name):
        try:
            return self._dim_names.index(dim_name)
        except ValueError:
            raise ValueError('"{}" is not a valid dimension name for StencilPoint'.format(dim_name)) from None

    def get_nonzero_dim_index(self, include_time=True):
        nz_dims = self.get_nonzero_dims(include_time=include_time)
        nz_dim_inds = tuple([self.get_dim_index(dim) for dim in nz_dims])
        return nz_dim_inds

    def get_dim_value(self, dim):
        return self.indices[self.get_dim_index(dim)]

    def change_dim(self, dim, new_value):
        if not isinstance(new_value, int):
            raise TypeError('new_value must be an int')

        dim_index = self.get_dim_index(dim)
        tmp_indices = list(self.indices)
        tmp_indices[dim_index] = new_value
        self.indices = tuple(tmp_indices)

    def swap_dim(self, new_dim):
        """
        Modify which spatial dimension the point lies along. The point must have only one nonzero spatial index, or this
        method will error
        :param new_dim: the dimension to change to. The index of the old nonzero dimension will be moved into this one.
        Must be one of the strings "x", "y", or "z"
        :return: none
        """
        nz_ind = self.get_nonzero_dim_index(include_time=False)
        if len(nz_ind) == 0:
            # No nonzero spatial indices, nothing to do
            return
        elif len(nz_ind) == 1:
            try:
                new_dim_idx = self._dim_names.index(new_dim)
            except ValueError:
                raise ValueError('new_dim must be one of "x", "y", or "z"') from None
            tmp_indices = list(self.indices)
            tmp_indices[new_dim_idx] = self.indices[nz_ind[0]]
            tmp_indices[nz_ind[0]] = 0
            self.indices = tuple(tmp_indices)
        else:
            raise RuntimeError('Cannot change the dimension of a point the spans >1 dimension')

    @staticmethod
    def expand_indices_to_dim(indices, dim=None):
        if dim is None:
            # Just need to make sure that there are four dimensions given. Append zeros as needed
            return indices + (0,) * (4 - len(indices))
        else:
            if len(indices) != 2:
                raise ValueError('If dim is specified, then indices must have only 2 elements')
            elif dim.lower() == 'x' or dim.lower == 'i':
                return indices[0], indices[1], 0, 0
            elif dim.lower() == 'y' or dim.lower == 'j':
                return indices[0], 0, indices[1], 0
            elif dim.lower() == 'z' or dim.lower == 'k':
                return indices[0], 0, 0, indices[1]
            else:
                raise ValueError('Value of dim "{}" not recognized. Allowed values are "x", "y", "z", "i", "j", or "k"'.format(dim))


class Stencil:
    """
    Stencil instances represent the points in space which are used to calculate a discrete derivative. For example, a
    centered first derivative (:math:`df/dx`) which is discretized as

    .. math::

        \frac{df}{dx} = \frac{f_{i+1} - f_{i-1}}{2 \Delta x}

    requires values from the i-1 and i+1 points, i.e. the point before and after the one at which we are calculating the
    derivative along the same dimension we are taking the derivative. Moreover, the i-1 point is subtracted from the i+1
    point. This would be represented by the following instance of Stencil:

        Stencil( (0,1), 1.0, (0,-1), -1.0, factored_coeff=1/(2*dx) )

    The first four arguments alternate indices and coefficients. Each index is given as a tuple, with the form
    (time_index, space_index). Both elements must be integers, and both are relative to the current point. The space
    index is assumed to be along the dimension that the derivative is with respect to. So (0, 1), 1.0 says that the
    point at the current time step (the 0 in the tuple) and 1 point further along in the differentiated dimension (the 1
    in the tuple, for this example, we're differentiating along the x dimension) should be multiplied by 1.0. Likewise,
    (0,-1), -1.0 means the point at the current timestep and one behind in the x dimension should be multiplied by -1.0.
    Finally, the factored_coeff keyword argument provides a way to describe a constant factor that all points should be
    multiplied by (i.e. the :math:`1/(2\Delta x)` here). If not given, this is set to 1.0.
    """

    #TODO: override __setattr__ to error if called and constant=True?
    def __init__(self, *args, factored_coeff=1.0, dim=None, constant=False):
        self._points = []
        self._constant = False  # must initialize, but set to the input value later after all the input points are added

        if len(args) % 2 != 0:
            raise ValueError('Stencil must be initialized with an even number of arguments, alternating indices and coefficients')

        for idx in range(0, len(args), 2):
            self.add_point_to_stencil(args[idx], args[idx+1]*factored_coeff, dim=dim)

        self._constant = constant  # safe to set now that the input points have been added

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError('Cannot add {} and {}'.format(self.__class__.__name__, other.__class__.__name__))

        new_stencil = self.duplicate()
        for point in other._points:
            new_stencil.add_point_to_stencil(point.indices, point.coefficient, action='add')

        return new_stencil

    def duplicate(self, new_time=None, new_dim=None, new_prefactor=1.0, keep_constant=False):

        point_args = [[pt.indices, pt.coefficient] for pt in self._points]
        point_args = itertools.chain.from_iterable(point_args)
        if keep_constant:
            constant = self._constant
        else:
            constant = False

        new_stencil = Stencil(*point_args, constant=constant, factored_coeff=new_prefactor)

        if new_time is not None:
            new_stencil.change_points_dim('t', new_time)
        if new_dim is not None:
            new_stencil.swap_dim(new_dim)

        return new_stencil

    def print_points(self):
        """
        Convenience method to print a human readable list of included points
        :return: none
        """
        slist = [str(pt) for pt in self._points]
        print('\n'.join(slist))

    def add_point_to_stencil(self, indices, coeff, action='replace', dim=None):
        """
        Adds a new coefficent to the internal list of coefficients.
        :param indices: A 2-element tuple consisting of (time_index, space_index).
        :param coeff: The coefficient the value of the point should be multiplied by in the discreet derivative.
        :param action: What to do if the point already has a coefficient. The default is 'replace', meaning that the
        existing coefficient is replaced. Alternately, 'add' will cause the new coefficient to be added to the existing
        one.
        :return: none
        """
        self._check_if_constant()

        # Does this point already exist?
        existing_point = self._find_point_if_defined(indices, dim=dim)
        if existing_point is None:
            self._points.append(StencilPoint(indices, coeff))
        elif action == 'replace':
            existing_point.coefficient = coeff
        elif action == 'add':
            existing_point.coefficient += coeff
        else:
            raise ValueError('No method implemented for action "{}"'.format(action))

    def change_points_dim(self, dim, new_value):
        if not isinstance(new_value, int):
            raise TypeError('new_time must be an int')

        dim_indices = [pt.get_dim_value(dim) for pt in self._points]
        check_same = [tidx == dim_indices[0] for tidx in dim_indices]
        if not all(check_same):
            raise RuntimeError('Cannot change dimension {} of points in stencil; values are not the same for all points'.format(dim))

        for point in self._points:
            point.change_dim(dim, new_value)

    def swap_dim(self, newdim):
        # Do all points in the stencil have the same nonzero dimension? And is there only one nonzero dimension?
        nz_dims = [pt.get_nonzero_dims(include_time=False) for pt in self._points if len(pt.get_nonzero_dims(include_time=False)) > 0]
        check_same = [d == nz_dims[0] and len(d) == 1 for d in nz_dims]
        if not all(check_same):
            raise RuntimeError('Cannot change spatial dimension of points in stencil; more than one spatial dimension is spanned')

        for point in self._points:
            point.swap_dim(newdim)

    def _find_point_if_defined(self, indices, dim=None):
        """
        Searches already stored points to find out if one at indices is already stored.
        :param indices: a 2-element tuple (time_index, space_index) or an instance of StencilPoint
        :return: the point instance if one exists, None otherwise.
        """

        if not isinstance(indices, tuple):
            test_point = StencilPoint(indices, dim=dim)
        elif isinstance(indices, StencilPoint):
            raise TypeError('indices must be an instance or tuple or StencilPoint')

        for point in self._points:
            if indices == point.indices:
                return point

        return None  # not necessary, since a function that does not return anything returns None, but make it explicit

    def construct_matrix(self, domain_size):
        # Construct the matrix row by row, adding each dimension's terms as needed.
        n_model_points = np.prod(domain_size)
        A = np.zeros((n_model_points, n_model_points))
        for row_idx, row in enumerate(A):
            self.apply_to_row(row, domain_size, row_idx)

        return A

    def apply_to_row(self, row, domain_size, curr_index):
        self._are_row_and_domain_consistent(row, domain_size)
        self._check_curr_index(curr_index, row)
        for point in self._points:
            if point.indices[0] == 0:
                # current timestep, needs to be added to the RHS vector, not LHS matrix
                pass
            elif point.indices[0] == 1:
                # next timestep, add to LHS matrix
                add_coefficient_to_row(row, point.coefficient, curr_index, domain_size,
                                       di=point.indices[1], dj=point.indices[2], dk=point.indices[3])
            else:
                raise NotImplementedError('Not configured to add points at timestep other than the current or next timestep (offending time index = {})'.format(point.indices[0]))


    @staticmethod
    def _are_row_and_domain_consistent(row, domain_size):
        if not isinstance(row, np.ndarray) or row.ndim != 1:
            raise TypeError('row must be numpy array')
        elif not isinstance(domain_size, tuple) or any([not isinstance(el, int) for el in domain_size]):
            raise TypeError('domain_size must be a tuple of ints')
        if row.size != np.prod(domain_size):
            raise ValueError('row.size must equal the product of the lengths of the dimensions given in domain_size')

        return True

    @staticmethod
    def _check_curr_index(curr_index, row):
        if not isinstance(curr_index, int):
            raise TypeError('curr_index must be an int')
        elif curr_index < 0 or curr_index >= row.size:
            raise IndexError('curr_index is out of range ({}-{}) for the given row'.format(0, row.size))

    def _check_if_constant(self):
        if self._constant:
            raise RuntimeError('Cannot modify a stencil set as constant. Make a copy with the duplicate() method and modify that.')

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


def domain_to_vector(domain):
    return domain.ravel(order='F')


def vector_to_domain(vector, domain_size):
    return vector.reshape(domain_size, order='F')

# Define common stencils for various ways of discretizing derivatives
# C.f. http://mathfaculty.fullerton.edu/mathews/n2003/differentiation/NumericalDiffProof.pdf
# tables 6.3 and 6.4
# These all assume that the factored coefficient will be set to D*dt/dx**n as required (any constants, i.e. how in a
# centered difference the 2 in the denominator of 2 * dx) are incorporated into the coefficients. Make all of these
# constant so that they can't be accidentally modified at run time.
# TODO: how to pass the grid spacing to these so that it is automatically raised to the appropriate power?
time_forward1_stencil = Stencil((1, 0), 1.0, (0, 0), -1.0, constant=True)
space_centered1_order2_stencil = Stencil((0, 1), 0.5, (0, -1), -0.5, constant=True)
space_centered1_order4_stencil = Stencil((0, 2), -1./12., (0, 1), 8./12., (0, -1), -8./12., (0, -2), 1./12., constant=True)
space_centered2_order2_stencil = Stencil((0, 1), 1.0, (0, 0), -2.0, (0, -1), 1.0, constant=True)
space_centered2_order4_stencil = Stencil((0, 2), -1./12., (0, 1), 16./12., (0, 0), -30./12., (0, -1), 16./12., (0, -2), -1./12., constant=True)