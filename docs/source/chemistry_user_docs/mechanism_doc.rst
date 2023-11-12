Chemistry in PECANS
===================

.. _ideal_chem_mech:

Idealized mechanisms
--------------------

**ideal_first_order**: This mechanism simply assumes that any chemical species are lost at the same first order rate.
That is, every species *A* follows the rate law:

.. math::
    \frac{\partial [A]}{\partial t} = -\frac{1}{\tau} [A]

:math:`\tau` is set in the configuration file using the ``mechanism_opts`` option. Specifically, it must include the
``lifetime_seconds`` key, e.g.::

    mechanism_opts = { lifetime_seconds = 10800 }

will set a 3 hour lifetime for any and all model species. The only other (optional) setting is ``species_name`` which 
will set the name of the model specie. If this is not given, the specie will just be called ``A``, and this is the 
specie name to use when setting emissions and initial conditions.

**ideal_two_phases_first_order**: This mechanism is similar to "ideal_first_order", except that it requires two lifetimes 
be specified. One lifetime will be in effect over the first part of the domain, the other over the second part. The options
for this mechanism are:

* ``first_lifetime_seconds``: the effective lifetime in seconds for the first part of the domain
* ``second_lifetime_seconds``: the effective lifetime in seconds for the second part of the domain
* ``first_stage_width``: this sets where the second lifetime takes effect. Specifically, PECANS will use the first lifetime 
  for all grid cells with :math:`x \leq x_\mathrm{emis} + w` where :math:`x_\mathrm{emis}` is the emissions center along the
  x-axis and :math:`w` is this option. (See below for a discussion of the required emissions settings.)
* ``species_name``: as with "ideal_first_order", the name of the sole specie in the model, defaulting to ``A``.

*Required emissions settings:* this mechanism **requires** that your configuration include the ``center_x`` option in the
``EMISSIONS`` section, even if you have emissions turned off. This contributes to setting where the second lifetime takes
effect. Note that for this mechanism you must use the "general" emissions configuration, where the settings are defined for
all species in the ``EMISSIONS`` section directly, i.e.::

    [EMISSIONS]
    emission_type = "gaussian"
    center_x = 10000
    width_x = 1000

.. _compiled_chem_mech:

Compiled mechanisms
-------------------

Currently, any non-ideal mechanism is selected by setting ``mechanism = "compiled"`` in the ``[CHEMISTRY]`` section. To switch
among different compiled mechanisms, you must run the ``build_pecans.py`` script to regenerate the Cython code for these mechanisms.
These mechanisms allow you to specify some species as "constant", meaning that their concentrations do not change over the course of 
the simulation. Each mechanism you wish to define must have a ``.spc`` and ``.eqn`` file in the :file:`pecans/Mechanisms` directory.
For details on these files, see :ref:`pecans_mech_format`.

Unlike the ideal mechanisms, compiled mechanisms do not use the ``mechanism_opts`` field, since there are not configuration options
for the mechanism itself. However, these mechanisms do accept constant parameters, which can be species concentrations or other 
required values that do not change over time. The easiest way to specify these values is if they should be the same for all model
grid cells. In that case, you can specify the value directly through the configuration. For example, here is the chemistry section
from the "one_box_nox_emis" example::

    [CHEMISTRY]
    do_chemistry = true
    mechanism = "compiled"
    initial_cond = "zero"
    const_params = {TEMP = 298.0, CAIR = 2e19, HO = 1e6, O3 = 8e11}

This is meant for use with the "nox" mechanism. Here, we specify that the temperature (``TEMP``) and number density of air 
(``CAIR``) should be 298 K and :math:`2 \times 10^{19}` molec. cm :math:`^{-3}`, respectively. These two parameters will be 
needed for all compiled mechanisms. In this case, we also want to fix the concentrations of two of our species, HO and O3.
Although the configuration file this example is taken from only has one grid cell in the model, this would also set these values for
all boxes in a multibox model.

.. note::
   Fixing the concentrations of one or more species can be useful if (1) the species is extremely short-lived and you are having
   trouble getting the resulting stiff system of differential equations to solve or (2) if you want to test the effect of a species
   being in a steady state.

If you need a parameter to be constant in time, but have different values in different model grid cells, you can put the string "file"
instead of a numeric value for a parameter. In this case you must also specify the ``const_param_input_file`` option to point to a file
that has these species. For example, the "100_box_nox_emis" example has::

    [DOMAIN]
    nx = 100
    ny = 0
    nz = 0

    [CHEMISTRY]
    do_chemistry = true
    mechanism = "compiled"
    initial_cond = "zero"
    const_params = {TEMP = 298.0, CAIR = 2e19, HO = "file", O3 = "file"}
    const_param_input_file = "examples/100_box_nox/100_box_nox_const.nc" 

In this case, the HO and O3 species will have their concentrations set by the arrays in the :file:`100_box_nox_const.nc` file. The 
contents of this file are::

    $ ncdump -h examples/100_box_nox/100_box_nox_const.nc
    netcdf \100_box_nox_const {
    dimensions:
        x = 100 ;
    variables:
        double HO(x) ;
            HO:units = "molec.cm^-3" ;
        double O3(x) ;
            O3:units = "molec.cm^-3" ;
    }

Notice that the file contains variables with the same names as the species we want to have time-constant, spatially-varying concentrations
and that these variables have the same shape as our domain, and they do *not* have a time dimension.

.. note::
   The path given for ``const_param_input_file`` in the config file can be absolute or relative. If relative, it must be relative to where
   you run PECANS from, *not* the location of the config file.

See also
--------

For the code used to set up the idealized solvers, see :ref:`ideal_chem_mech_code`
