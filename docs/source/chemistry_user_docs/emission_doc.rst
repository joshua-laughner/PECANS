.. _emission_doc:

Emissions options
=================

Emissions in PECANS can be defined once for all species (which is usually only used for the idealized mechanisms)
or on a per-species basis. At the moment, emissions can only be given predefined shapes, there is no option to input
them through a netCDF file. When defining the emissions for all species, the options for the chosen emission type
are given directly in the emissions section::

  [EMISSIONS]
  do_emissions = true
  emission_type = "point"
  center_x = 0.0
  center_y = 0.0
  center_z = 0.0
  total = 1e24


To define per-species emissions, the options for each specie's emissions must go in ``[[EMISSIONS.species]]`` sections
like so::

  [EMISSIONS]
  do_emission = true

  [[EMISSIONS.species]]
  emission_specie = "NO"
  emission_type = "point"
  ...

  [[EMISSIONS.species]]
  emissions_specie = "NO2"
  emission_type = "gaussian"
  ...

(Note that the ``...`` indicates options omitted for brevity.) Any specie without a corresponding section will have zero
emissions. 

The following sections describe the required options for each emission type.

gaussian
--------

Emissions are distributed in a 1D, 2D, or 3D Gaussian. Its shape are center are determined by the following options:

- *center_x*, *center_y* - the center coordinates (in meters) of the Gaussian. If 1D, *center_y* may be omitted.
- *width_x*, *width_y* - the width (:math:`\sigma`) of the Gaussian in each dimension. The *y*
  width may be omitted according to the same rules as the centers.
- *total* - the total emissions, in molecules per second. Note that this is different from the Gaussian used for
  initial conditions, because this sets the total area under the Gaussian, not its maximum value.

In 1D, the emissions at each point would be given by:

.. math::
   E(x) = \frac{E_{\mathrm{tot}} \Delta x}{\Delta x \Delta y \sqrt{2\pi\sigma_x^2}} \exp\left( -\frac{x - \mu_x}{2\sigma^2}\right)

where :math:`\Delta x` and :math:`\Delta y` are the box ``dx`` and ``dy`` values, :math:`E_{\mathrm{tot}}` is the total
emissions specified by the *total* emissions option, :math:`\mu_x` is the *center_x* and :math:`\sigma_x` the *width_x*.
The total emissions are divided by :math:`\Delta x \Delta y` to normalize per unit area, but must be multiplied by
:math:`\Delta x` to account for the discretization of the Gaussian onto the model grid (i.e. the normalized Gaussian is
the emissions fraction per unit length, so it must be multiplied by the width of the box).

Note that Gaussian emissions will always be placed at the surface.


point
-----
Emissions are all placed in a single grid cell. Its quantity and location are determined by the following options:

- *center_x*, *center_y*, *center_z* - the center coordinates (in meters) where the emissions should be placed.
  Whichever grid cell this falls in received the emissions. If 1D, *center_y* and *center_z* may be omitted; if 2D,
  *center_z* may be omitted.
- *total* - the total emissions, in molecules per second.

Note that unlike Gaussian emissions, point emissions can be placed in any box (surface or not).

constant
--------
The same emissions are set for all grid cells. Its quantity is determined by the following options:

- *rate* - the emission rate for *each* grid cell, in molecules per second.
- *surface_only* - by default, in a 3D model, these emissions will only be placed in the surface boxes.
  To override that and place them in all boxes, set this option to ``false``.