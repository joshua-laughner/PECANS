PECANS Configuration Manual
===========================

The config file
---------------

The configuration file for PECANS is a `fairly standard format <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`_
for Linux-y configuration files. Sections are demarcated by all-caps text inside brackets; options within each section
are given starting at the beginning of the line, followed by an equals sign, then the option value.

The values will get parsed into Python literals, following certain rules:

* Simple values:

    + **booleans**: written as either ``True`` or ``False``, with that exact capitalization
    + **integers**: numbers without a decimal point, may start with a ``+``, ``-``, or number. E.g., ``1``, ``-20``, and
      ``+300`` are all valid integers
    + **floating point numbers** a.k.a. *real numbers*: numbers with a decimal point, may start with or without ``+`` or ``-``.
      E.g., ``1.``, ``1.0``, and ``+1.0`` are all the same value.
    + **strings**: anything not matching these rules (and not a compound value) is kept as a string

* Compound values:

    + **tuples**: if a value needs to contain multiple values, then separate each value with a comma. This will be stored
      in Python as a tuple (basically an immutable list).
    + **dictionaries**: basically a way to name sub options instead, these need to be written as ``key1: value1, key2: value2``
      i.e. the name of each sub option, followed by a colon, then the value. The value must be one of the simple value types.

Most of these will be used somewhere in the defaults, so peruse the default configuration for examples.

Output
------

* ``output_frequency`` - how often (in seconds) model output should be written. This does not need to be a multiple of the
  time step, the model will write output if the time since the last output file is greater than this number. This does
  mean that if this is *not* a multiple of the time step that the spacing between files may be inconsistent


Domain
------

Options affecting domain size
*****************************

* ``nx`` - number of model boxes in the *x* dimension. Must be a scalar integer > 0.
* ``ny`` - number of model boxes in the *y* dimension. Must be a scalar integer >= 0. Unlike ``nx``, this may be 0, which
  means that the model will be 1D.
* ``nz`` - number of model boxes in the *z* dimension. Must be a scalar integer >= 0. If ``ny`` is 0, this must also be 0;
  if this is 0 and ``ny`` is not 0, then the model will be 2D.

* ``dx`` - size of the boxes in the *x* dimension in meters. Must be a scalar integer > 0.
* ``dy`` - as ``dx`` but for the *y* dimension. Must be >0 even for 1D models.
* ``dz`` - as ``dx`` but for the *z* dimension. Must be >0 even for 1D or 2D models.

Options affecting model time
****************************

* ``dt`` - model timestep in seconds. The upper limit on this will be set by the stability of the chemistry, transport,
  etc.
* ``run_time`` - how long (in seconds) the model should run for.


Transport
---------

* ``do_transport`` - a boolean (``True`` or ``False``) that turns transport calculation on or off.
* ``scheme`` - what solver to use for transport. Options are:

    + **implicit2** or **backwards_euler_2** - a second order implicit solver. Theoretically unconditionally stable
      (meaning that your ``dt`` may be large) but potentially less accurate and more computationally expensive than other
      methods
    + **crank_nicholson** - a higher-order method that is not unconditionally stable, but may be more accurate than a
      fully implicit method.

* ``wind_type`` - determines how winds are calculated for each time step. Options are:

    + **fixed** - the winds and turbulent diffusion constants in the x, y, and z directions are specified in the
      configuration file (see ``wind_speeds`` and ``diffusion_coeffs``).

* ``wind_speeds`` - a dictionary of the form ``x: U, y: V, z: W`` where *U*, *V*, and *W* are the wind speeds in meters per
  second along the *x*, *y*, and *z* directions, respectively. Can be omitted is ``wind_type`` is not "fixed".
* ``diffusions_coeffs`` - a dictionary of the form ``x: Dx, y: Dy, z: Dz`` where *Dx*, *Dy*, and *Dz* are the diffusion
  coefficients in meters squared per second in the *x*, *y*, and *z* directions, respectively. (100 is a good default
  value for a moderately turbulent atmosphere). Can be omitted if ``wind_type`` is not "fixed".


Chemistry
---------

* ``do_chemistry`` - a boolean (``True`` or ``False``) that turns chemistry on or off.
* ``mechanism`` - determines which chemical mechanism the model uses. There are both idealized mechanisms and explicit
  mechanisms. Idealized options are:

    + **ideal_first_order** - all chemical species are removed with a characteristic lifetime

* ``mechanism_opts`` - additional options to be set for each mechanism, as a dictionary. For details, see
  :ref:`ideal_chem_mech`.
* ``initial_cond`` - determines how the initial chemical concentrations are set. Options are:
    + **gaussian** - sets the initial conditions for all species as a Gaussian with a given center, width, and height.
      These values are set by the ``initial_cond_opts`` line.
* ``initial_cond_opts`` - additional options required by whatever initial conditions are selected. See
  :ref:`initial_cond_doc` for information on the specific options required for each initial condition.


Emissions
---------

* ``do_emissions`` - a boolean (``True`` or ``False``) that turns emissions on or off.
* ``emission_type`` - determines how the emissions are distributed throughout the domain. Options are:
    + **gaussian** - emissions are distributed in a 1D, 2D, or 3D Gaussian as appropriate.