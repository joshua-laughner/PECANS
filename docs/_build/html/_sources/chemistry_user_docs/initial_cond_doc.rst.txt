.. _initial_cond_doc:

Initial conditions
==================

Here is a list of valid options for the ``initial_cond`` option in the configuration file:

* **gaussian** - Initialize all chemical species to be a 1D, 2D, or 3D Gaussian as appropriate for the model domain.
  Requires the following options for ``initial_cond_opts``:

    + *center_x*, *center_y*, *center_z* - the center point of the Gaussian in the *x*, *y*, and *z* dimensions, given
      in meters. If the model is 1D, *center_y* and *center_z* may be omitted; likewise *center_z* may be omitted if
      the model is 2D.
    + *width_x*, *width_y*, *width_z* - the width of the Gaussian (its :math:`\sigma`) in meters for the *x*, *y*, and *z*
      dimensions, respectively. Like the center coordinates, the y and z may be omitted if the model doesn't include that
      dimension.
    + *max* - the maximum concentration, in molecules per cubic centimeter.

  In 1D that means that the initial concentrations of each species, :math:`A` is given by:

  .. math::
    [A](x) = m_A \exp\left(-\frac{x - \mu_x}{2\sigma_x^2}\right)

  where :math:`m_a` is the maximum given by the *max* option, :math:`\mu_x` is *center_x* and :math:`\sigma_x` is
  *width_x*