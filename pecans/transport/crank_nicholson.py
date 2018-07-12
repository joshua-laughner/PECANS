from . import transport_utils as tutils

def construct_transport_matrix_with_stencil(dt, dx=None, dy=None, dz=None, u_x=None, u_y=None, u_z=None, D_x=None, D_y=None, D_z=None, domain_size=None, boundary_conditions=None):
    # Returning these values from the check function ensures that they are converted to the proper types (i.e. ints to
    # floats) and default values are assigned.
    dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions, n_model_dims = \
        tutils.check_transport_inputs(dt, dx, dy, dz, u_x, u_y, u_z, D_x, D_y, D_z, domain_size, boundary_conditions)


    # Construct the combined quantities used in the matrix
    r_x = D_x * dt / (dx**2)
    C_rx = u_x * dt / dx  # the absent factor of 2 in the denominator is handles by the stencil

    # Construct the Stencil instances that represent the time derivative, advection derivative, and diffusion derivative
    time_stencil = tutils.time_forward1_stencil.duplicate()
    advection_stencil = tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_rx*0.5, new_time=1) \
                        + tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_rx*0.5, new_time=0)
    diffusion_stencil = tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_x*0.5, new_time=1) \
                        + tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_x*0.5, new_time=0)

    if n_model_dims >= 2:
        r_y = D_y * dt / (dy**2)
        C_ry = u_y * dt / dy
        advection_stencil += tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_ry*0.5, new_time=1, new_dim='y') \
                             + tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_ry*0.5, new_time=0, new_dim='y')
        diffusion_stencil += tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_y*0.5, new_time=1, new_dim='y') \
                             + tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_y*0.5, new_time=0, new_dim='y')

    if n_model_dims >= 3:
        r_z = D_z * dt / (dz**2)
        C_rz = u_z * dt / dz
        advection_stencil += tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_rz*0.5, new_time=1, new_dim='z') \
                             + tutils.space_centered1_order2_stencil.duplicate(new_prefactor=C_rz*0.5, new_time=0, new_dim='z')
        diffusion_stencil += tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_z*0.8, new_time=1, new_dim='z') \
                             + tutils.space_centered2_order2_stencil.duplicate(new_prefactor=-r_z*0.8, new_time=0, new_dim='z')

    total_stencil = time_stencil + advection_stencil + diffusion_stencil

    return total_stencil.construct_matrix(domain_size)