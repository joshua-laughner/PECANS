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

    mechanism_opts = lifetime_seconds: 10800

will set a 3 hour lifetime for any and all model species.

For the code used to set up these solvers, see :ref:`ideal_chem_mech_code`
