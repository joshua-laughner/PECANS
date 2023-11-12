.. _initial_cond_doc:

Initial conditions
==================

The ``initial_cond`` option in a PECANS configuration file can be either the string "zero" or a dictionary
of dictionaries specifying initial conditions for each specie in the model. Specifying the string "zero" will
make all species start with a concentration of 0.0. Otherwise, you can specify each specie to have an initial 
condition of "gaussian", "point", or "flat":

* "gaussian" will create an initial condition with a distribution that follows a 1D, 2D, or 3D Gaussian shape
  (depending on the model dimensions).
* "point" will create an initial condition with the concentration all in a single grid cell.
* "flat" will create an initial condition with the same concentration in all grid cells.

An example using each of these options for a 1D model is::

  [CHEMISTRY]
  do_chemistry = true
  mechanism = "compiled"
  const_params = {TEMP = 298.0, CAIR = 2e19, HO = 1e6}

  [[CHEMISTRY.initial_cond]]
  specie = "NO"
  initial_type = "gaussian"
  max_concentration = 1e10
  center_x = 50000
  width_x = 5000

  [[CHEMISTRY.initial_cond]]
  specie = "NO2"
  initial_type = "point"
  center_x = 50000
  concentration = 4e10

  [[CHEMISTRY.initial_cond]]
  specie = "O3"
  initial_type = "flat"
  concentration = 1e12

Notice that we repeat the ``[CHEMISTRY.initial_cond]]`` header for each specie that needs an initial condition.
The double square bracket notation is how we define a list at the top level of a TOML file. Any specie not 
listed will have an initial concentration of zero.

.. note::
   For the ideal mechanisms, the sole specie is named "A" by default, so setting up initial conditions for that should
   use ``specie = "A"`` unless you change the specie name.

The following sections describe each of the initial condition types in detail.

Gaussian
--------
Initialize the chemical specie to be a 1D, 2D, or 3D Gaussian as appropriate for the model domain.
Requires the following options:

- *center_x*, *center_y*, *center_z* - the center point of the Gaussian in the *x*, *y*, and *z* dimensions, given
  in meters. If the model is 1D, *center_y* and *center_z* may be omitted; likewise *center_z* may be omitted if
  the model is 2D.
- *width_x*, *width_y*, *width_z* - the width of the Gaussian (its :math:`\sigma`) in meters for the *x*, *y*, and *z*
  dimensions, respectively. Like the center coordinates, the y and z may be omitted if the model doesn't include that
  dimension.
- *max_concentration* - the maximum concentration, in molecules per cubic centimeter.

In 1D that means that the initial concentrations of each species, :math:`A` is given by:

.. math::
   [A](x) = m_A \exp\left(-\frac{x - \mu_x}{2\sigma_x^2}\right)

where :math:`m_a` is the maximum given by the *max_concentration* option, :math:`\mu_x` is *center_x* and :math:`\sigma_x` is
*width_x*.


Point
-----
Initializes the chemical specie to have a given concentration in a single grid cell. Requires the following options:

- *center_x*, *center_y*, *center_z* - defines which grid cell to place the initial concentration in. These give the center in meters
  from the origin. Whichever grid cell contains this coordinates will have the initial concentration. If the model is 1D, *center_y* 
  and *center_z* may be omitted; likewise *center_z* may be omitted if the model is 2D.
- *concentration* - defines the initial concentration to give the cell in molec. cm :math:`^{-3}`.

.. warning::
   Because the point initial condition creates extremely sharp gradients in the model concentrations, you may see unphysical osciallations
   in the output concentrations. Generally a gaussian initial condition is better if you need to start with concentrations in a specific location.

Flat
----
Initializes the chemical specie to have a given concentration in all grid cells. Requires the following option:

- *concentration* - defines the initial concentration to give all grid cells in molec. cm :math:`^{-3}`.