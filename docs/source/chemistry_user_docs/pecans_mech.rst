.. _pecans_mech_format:

PECANS Mechanism Format
=======================

PECANS has a format for specifying compiled mechanisms that is modeled on the Kinetic Pre-processor (KPP) format,
but somewhat simplified. It consists of three file types: species files (``.spc``), equation/reaction files (``.eqn``),
and the rate expression files (``.rate``).

Species files
-------------

For PECANS, a species file is just a file listing the species in the mechanism, with one specie per line. Our simple
NOx mechanism has the species file::

    NO
    NO2
    O3
    HO
    HNO3

The main purpose of the species file is to help catch typos, as it is an error if an equation file uses a specie not
listed in the species file.

Equation file
-------------

The equation (a.k.a. reaction) file specifies the chemical reactions to simulate. An example from the simple NOx example
mechanism is::

    NO2 = NO + O3 : 0.01
    NO + O3 = NO2 : ARR2( 1.40e-12, 1310.0, TEMP )
    HO + NO2 = HNO3 : TROE( 1.49e-30 , 1.8, 2.58e-11 , 0.0, TEMP, CAIR)

To break this down, each line has three parts: reactants, products, and the rate expression. Reactants and products are separated
by the ``=`` and the rate expression follows the ``:``. If a product or reactant needs a coefficient, that coefficient must go in
front of the specie separated by a space, e.g.::

    H2O2 = 2 HO

The rate expression must be a valid Python expression. In the above example, we see two cases. The simplest is the ``NO2 = NO + O3``
reaction, which just has a number. This will produce a reaction with a constant rate.

.. note::
   The ``NO2 = NO + O3`` reaction has its rate expression as a constant because this is a photolysis reaction and PECANS does not 
   have the ability to calculate photolysis rates yet. In principle, these could be parameterized as rate expressions that depend on
   solar zenith angle or other proxy variables.

The second and third equations use functions defined in the rates file (see next section). Usually,
kinetic rate expressions follow one of a few forms with different constants, and depend on temperature and/or the number density of air.
For ``NO + O3 = NO2``, we use the ``ARR2`` function, which is defined in the rates file as::

    cdef ARR2(double A0, double B0, double TEMP):
        return A0 * exp(-B0 / TEMP)

We can see that it takes two constants (``A0`` and ``B0``) and temperature. This matches the call in the equation file, where two constants
are given explicitly (``A0 = 1.40e-12``, ``B0 = 1310.0``) and temperature is given as the variable ``TEMP``.

.. note::
   PECANS uses ``TEMP`` for temperature and ``CAIR`` for number density of air in both the mechanisms and configuration. This is hard-coded
   to allow it to check that these parameters are included in the configuration if needed.

We will discuss how rates are defined :ref:`below <mech_rates_file>`.

Some mechanisms will want to use additional parameters in the rate constants. For example, the ``nox_voc`` mechanism uses ``alpha`` 
to represent the branching ratio of some RO2 + NO reactions. At present, these constants must be specified when the mechanism is
built using the ``--params`` keyword for :file:`build_pecans.py`, for example::

    $ python build_pecans.py --params alpha=0.05

Multiple parameters can be specified in sequence as::

    $ python build_pecans.py --params alpha=0.05 beta=10.0

To pass a mechanism name as a positional argument along with the ``--params`` arguments, either pass the mechanism name first or
use ``--`` to separate the position argument from the ``--params`` values::

    # Either will work
    $ python build_pecans.py nox --params alpha=0.5
    $ python build_pecans.py --params alpha=0.5 -- nox

The main limitation of the rate definition in the equation file is that each definition can only use one function defined in the 
:ref:`rates files <mech_rates_file>`.

.. _mech_rates_file:

Rates files
-----------

These files go in the :file:`pecans/Rates` folder and have the extension ``.rate``. Unlike the species and equations files, the rates
files are not associated with any particular mechanism. Rather, all ``.rate`` files in the Rates folder can be used by any mechanism.
These files simply define Cython ``cdef`` functions that calculate a kinetic rate given temperature, number density of air, and any number
of constants. See the `Cython language basics <https://cython.readthedocs.io/en/latest/src/userguide/language_basics.html>`_ for an 
introduction to writing Cython code.

See also
--------

Source of kinetic rates include the `JPL Data Evaluation <https://jpldataeval.jpl.nasa.gov/>`_ and the 
`Master Chemical Mechanism <https://mcm.york.ac.uk/MCM>`_.