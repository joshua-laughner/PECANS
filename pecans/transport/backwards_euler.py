import numpy as np

from . import transport_utils as tutils

import pdb

def construct_transport_equation_2nd_order_centered_space(dt, dx=None, dy=None, dz=None, u_x=None, u_y=None, u_z=None, D_x=None, D_y=None, D_z=None, domain_size=None, boundary_conditions=None):
    """
    Contructs the matrix of coefficients necessary to solve the advection-diffusion equation using
    a backwards euler method.
    :param dt: model time step in seconds. Must be a scalar number.
    :param dx: model spacing in the x direction in meters. Must be a scalar number (unequal spacing not implemented).
    :param dy: model spacing in the y direction in meters. Must be a scalar number (unequal spacing not implemented), or
    omitted for a 1D model.
    :param dz: model spacing in the z direction in meters. Must be a scalar number (unequal spacing not implemented), or
    omitted for a 1D or 2D model.
    :param u_x: The wind in the x-direction (meters/second), positive is to the east. This must either be a scalar number
    or a numpy array. If a scalar value is given, then the wind is assumed to be the same at all grid points
    in the model. In this case, the domain_size parameter must also be given. If given as a numpy array, then
    it must be 1D if only u_x and D_x are given, 2D if u_y and D_y are given, and 3D if u_z and D_z are given.
    To avoid unintended effects, all given u and D inputs must be either scalars or numpy arrays; you cannot
    mix the two. If given as numpy arrays, all must be the same size.
    :param u_y: The wind in the y-direction (meters/second), positive is to the north. If omitted, the model must be 1D.
    :param u_z: The wind in the z-direction (meters/second), positive is to up. If omitted, the model must be 1D or 2D.
    :param D_x: The diffusion constant in the x direction (meters squared per second). Must be positive.
    :param D_y: The diffusion constant in the y direction (meters squared per second). Must be positive. If omitted, the
    model must be 1D.
    :param D_z: The diffusion constant in the z direction (meters squared per second). Must be positive. If omitted, the
    model must be 1D or 2D.
    :return: a numpy matrix (A) and vector (b) to be solved as in Ax = b, where x is the concentrations at the next time
    step

    Details:
        In a backward Euler (a.k.a. implicit) approach, the spatial derivatives are evaluated at the next time step,
    rather than the current one. So, in the 1D advection diffusion equation:

    .. math::

        \frac{\partial c}{\partial t} + u_x \frac{\partial c}{\partial x} = D_x \frac{\partial^2 c}{\partial x^2}

    both derivatives with respect to x are evaluated at timestep n+1. In the following, :math:`i` subscripts are spatial
    indices and :math:`n` superscripts are time indices.

    .. math::
        \frac{c_i^{n+1} - c_i^n}{\Delta t} + u_x \frac{c_{i+1}^{n+1} - c_{i-1}^{n+1}}{2\Delta x} = D_x \frac{c_{i+1}^{n+1} - 2c_i^{n+1} + c_{i-1}^{n+1}}{\Delta x^2}

    This uses a centered difference scheme for both spatial derivatives. If we let :math:`r_x = D_x \Delta t / \Delta x^2`
    and :math:`C_{r,x} = u_x \Delta t / (2\Delta x)`, then we can rearrange to place all terms at timestep :math:`n+1`
    on the LHS and terms at timestep :math`n` on the RHS:

    .. math::

        (C_{r,x} - r)c_{i+1}^{n+1} + (1 + 2r)c_i^{n+1} + (-C_{r,x} - r)c_{i-1}^{n+1} = c_i^n
        :label:advdiff1D

    This can be transformed into a matrix equation :math:`Ax = b` where :math:`x = c^{n+1}` and :math:`b = c^n`.
    :math:`A` then is a matrix with coefficients of :math:`1 + 2r` along the main diagonal, :math:`C_{r,x}` - r on the
    upper diagonal and :math:`-C_{r,x} - r` along the lower diagonal.

    In 2D, we need to add :math:`u_y \partial c/\partial y` and :math:`D_y \partial^2 c/\partial y^2` and
    likewise for :math:`z` in 3D. Each of these results in an equation similar to :eq:`advdiff1D`, except the spatial
    indices are modified over the :math:`y` or :math`z` dimension. Adding a second dimension results in additional
    diagonals being filled in the matrix as coefficients, making it a pentadiagonal matrix. It will also add an extra
    factor of 2r to the main diagonal, because the :math:`c_i^{n+1} \equiv c_{ijk}^{n+1}` point shows up in the second
    derivative for diffusion in all three dimensions.

    In order to write a 2D or 3D domain in this matrix equation form, we need to "flatten" the domain into a vector.
    Essentially, each grid point is ordered according to its linear index. Concretely, if we had a :math:`3 \times 3` 2D
    domain, with points :math:`c_{11}` to :math:`c_{33}`. Both the :math:`c^n` and :math:`c^{n+1}` vectors would be:

    .. math::

        \begin{bmatrix} c_{11} \\ c_{21} \\ c_{31} \\ c_{12} \\ c_{22} \\ c_{32} \\ c_{13} \\ c_{23} \\ c_{33} \end{bmatrix}

    The exact order is unimportant, so long as the mapping from the 2- or 3-D space to the vector is consistent.
    """

    # Returning these values from the check function ensures that they are converted to the proper types (i.e. ints to
    # floats) and default values are assigned.
    dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions, n_model_dims = \
        tutils.check_transport_inputs(dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions)

    # Construct the combined quantities used in the matrix
    r_x = D_x * dt / (dx**2)
    C_rx = u_x * dt / (2*dx)

    if n_model_dims >= 2:
        r_y = D_y * dt / (dy**2)
        C_ry = u_y * dt / (2*dy)
    else:
        r_y = 0.0
        C_ry = 0.0

    if n_model_dims >= 3:
        r_z = D_z * dt / (dz**2)
        C_rz = u_z * dt / (2*dz)
    else:
        r_z = 0.0
        C_rz = 0.0

    # Construct the matrix row by row, adding each dimension's terms as needed.
    n_model_points = np.prod(domain_size)
    A = np.zeros((n_model_points, n_model_points))
    for i_row, row in enumerate(A):
        # 1D coefficients are always required. For consistency and simplicity, we'll use numpy's ravel and
        # ravel_multi_index methods
        tutils.add_coefficient_to_row(row, 1 + 2*r_x, i_row, domain_size)
        tutils.add_coefficient_to_row(row, C_rx - r_x, i_row, domain_size, di=1)
        tutils.add_coefficient_to_row(row, -C_rx - r_x, i_row, domain_size, di=-1)

        if n_model_dims >= 2:
            tutils.add_coefficient_to_row(row, 2*r_y, i_row, domain_size)
            tutils.add_coefficient_to_row(row, C_ry - r_y, i_row, domain_size, dj=1)
            tutils.add_coefficient_to_row(row, -C_ry - r_y, i_row, domain_size, dj=-1)

        if n_model_dims >= 3:
            tutils.add_coefficient_to_row(row, 2*r_z, i_row, domain_size)
            tutils.add_coefficient_to_row(row, C_rz - r_z, i_row, domain_size, dk=1)
            tutils.add_coefficient_to_row(row, -C_rz - r_z, i_row, domain_size, dk=-1)

    return A


def construct_transport_matrix_with_stencil(dt, dx=None, dy=None, dz=None, u_x=None, u_y=None, u_z=None, D_x=None, D_y=None, D_z=None, domain_size=None, boundary_conditions=None):
    # Returning these values from the check function ensures that they are converted to the proper types (i.e. ints to
    # floats) and default values are assigned.
    dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions, n_model_dims = \
        tutils.check_transport_inputs(dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions)

    # TODO: this breaks if given an array for the u's or D's. It's going to have to create an array of stencils or
    # something to handle that case. Or separate out the prefactor more carefully in the stencil implementation

    # Construct the combined quantities used in the matrix
    r_x = D_x * dt / (dx**2)
    C_rx = u_x * dt / dx  # the absent factor of 2 in the denominator is handles by the stencil

    # Construct the Stencil instances that represent the time derivative, advection derivative, and diffusion derivative
    time_stencil = tutils.time_forward1_stencil.duplicate()
    advection_stencil = tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_rx, new_time=1)
    diffusion_stencil = tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_x, new_time=1)

    if n_model_dims >= 2:
        r_y = D_y * dt / (dy**2)
        C_ry = u_y * dt / dy
        advection_stencil += tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_ry, new_time=1, new_dim='y')
        diffusion_stencil += tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_y, new_time=1, new_dim='y')

    if n_model_dims >= 3:
        r_z = D_z * dt / (dz**2)
        C_rz = u_z * dt / dz
        advection_stencil += tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_rz, new_time=1, new_dim='z')
        diffusion_stencil += tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_z, new_time=1, new_dim='z')

    total_stencil = time_stencil + advection_stencil + diffusion_stencil

    return total_stencil.construct_matrix(domain_size)
