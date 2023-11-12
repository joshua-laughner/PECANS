.. _config_user_guide:


PECANS Configuration Manual
===========================

The config file
---------------

Starting with v0.2.0, the configuration file for PECANS is a `TOML file <https://toml.io/en/>`_. We switched away from
the previous INI file type because TOML automatically handles the nested lists/dictionaries that we had to manually
parse before. Most usage should be pretty evident from the examples. For details on TOML syntax, see
https://toml.io/en/v1.0.0.

.. note::
   If you have existing v0.1 PECANS configuration files you want to convert, the main differences in syntax are:
     1. boolean values must be lowercase in TOML (e.g. ``true``), whereas before they had to be capitalized (``True``)
     2. string *values* must be quoted; strings used as keys can still be unquoted.
     3. inline dictionaries (note that TOML calls dictionaries "tables") use the syntax ``{ key1 = value1, key2 = value2, ...}``,
        whereas before we used ``key1: value1, key2: value2, ...``.

   This only covers the difference in INI-style syntax vs. TOML-style syntax. There were also significant changes to how the
   :ref:`initial conditions <initial_cond_doc>` and :ref:`emissions <emission_doc>` are specified, so read those sections of the documentation carefully if you have old configuration files you need to update.

Output
------

* ``output_frequency`` - how often (in seconds) model output should be written. This does not need to be a multiple of the
  time step, the model will write output if the time since the last output file is greater than this number. This does
  mean that if this is *not* a multiple of the time step that the spacing between files may be inconsistent
* ``output_path`` - where to save the output files. If not given, defaults to the current directory.

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
    + **ideal_two_phases_first_order** - all chemical species are removed with a characteristic lifetime that changes at a selected point in the domain.
    + **compiled** - chemistry uses the mechanism compiled by running ``build_pecans.py``.

* ``mechanism_opts`` - additional options to be set for each mechanism, as a dictionary. For details, see
  :ref:`ideal_chem_mech`.
* ``initial_cond`` - determines how the initial chemical concentrations are set. See :ref:`initial_cond_doc` for details.
* ``const_params`` - a dictionary specifying chemical species or other parameters that should not vary with time
  in a compiled mechanism. See the :ref:`compiled_chem_mech` section for details.
* ``const_param_input_file`` - a path to a netCDF file containing concentrations of chemical species or values of other
  parameters that will be constant in time. See :ref:`compiled_chem_mech` for details.

Emissions
---------

* ``do_emissions`` - a boolean (``True`` or ``False``) that turns emissions on or off.

Other options in this section differ depending on whether you specify per-specie emissions or identical emissions for all 
species. See :ref:`emission_doc` for details.