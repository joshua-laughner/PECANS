.. _emission_doc:

Emissions options
=================

* **gaussian** - Emissions are distributed in a 1D, 2D, or 3D Gaussian. Its shape are center are determined by the
  values given in ``emission_opts``:

    + *center_x*, *center_y*, *center_z* - the center coordinates (in meters) of the Gaussian. If 1D, *center_y* and
      *center_z* may be omitted; if 2D, *center_z* may be omitted.
    + *width_x*, *width_y*, *width_z* - the width (:math:`\sigma`) of the Gaussian in each dimension. The *y* and *z*
      widths may be omitted according to the same rules as the centers.
    + *total* - the total emissions, in molecules per second. Note that this is different from the Gaussian used for
      initial conditions, because this sets the total area under the Gaussian, not its maximum value.

  In 1D, the emissions at each point would be given by:

  .. math::
    E(x) = \frac{E_{\mathrm{tot}} \Delta x}{\Delta x \Delta y \sqrt{2\pi\sigma_x^2}} \exp\left( -\frac{x - \mu_x}{2\sigma^2}\right)

  where :math:`\Delta x` and :math:`\Delta y` are the box ``dx`` and ``dy`` values, :math:`E_{\mathrm{tot}}` is the total
  emissions specified by the *total* emissions option, :math:`\mu_x` is the *center_x* and :math:`\sigma_x` the *width_x*.
  The total emissions are divided by :math:`\Delta x \Delta y` to normalize per unit area, but must be multiplied by
  :math:`\Delta x` to account for the discretization of the Gaussian onto the model grid (i.e. the normalized Gaussian is
  the emissions fraction per unit length, so it must be multiplied by the width of the box).

* **point** - Emissions are all placed in a single grid cell. Its quantity and location are determined by the values
  given in ``emission_opts``:

    + *center_x*, *center_y*, *center_z* - the center coordinates (in meters) where the emissions should be placed.
      Whichever grid cell this falls in received the emissions. If 1D, *center_y* and *center_z* may be omitted; if 2D,
      *center_z* may be omitted.
    + *total* - the total emissions, in molecules per second.