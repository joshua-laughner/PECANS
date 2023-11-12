PECANS Ensemble Runners
=======================

PECANS includes the capability to run an ensemble of models, all based on an initial configuration,
but with changes made for each ensemble member. This is done through the :mod:`pecans.ensembles.api`
module, documented below. You will need to import this API and write a script using it; there is no
command line interface for it.

.. _ensemble_tutorial:

Tutorial
--------

To begin import the API in your script and create and instance of the :class:`~pecans.ensembles.api.EnsembleRunner`::

    from pecans.ensembles.api import EnsembleRunner

    runner = EnsembleRunner(
        base_config_file = 'examples/one_box_ideal/one_box_ideal.toml',
        ensemble_mode='iterations',
        root_output_dir='test_output'
    )

This has created a runner that will use the "one_box_ideal" example as the baseline configuration and which will
output its data under the :file:`test_output` directory. (Note, this assumes we will run this script from the root
of the PECANS repo, so these paths are relative to that.) We'll come back to the ensemble mode in a bit.

Next we need to define what options should change for each member of the ensemble. Let's say that we wanted to
test the effect of different lifetimes and initial concentration on the output. To do that, we would do::

    runner.add_ens_var_by_string('CHEMISTRY/mechanism_opts/lifetime_seconds', [3600, 7200, 10800])
    runner.add_ens_var_by_string('CHEMISTRY/initial_cond/0/concentration', [2, 20, 200])

The first argument to this method is the path to the option we want to vary. The ``CHEMISTRY`` section of our
config file is::

    [CHEMISTRY]
    do_chemistry = true
    mechanism = "ideal_first_order"
    mechanism_opts = {lifetime_seconds = 3600}

    [[CHEMISTRY.initial_cond]]
    specie = "A"
    initial_type = "point"
    center_x = 500
    concentration = 1

so the first option, ``'CHEMISTRY/mechanism_opts/lifetime_seconds'`` contains the keys for each of the dictionaries we need
to access separated by slashes: "CHEMISTRY" in the top dictionary, "mechanism_opts" in the ``CHEMISTRY`` dict, and 
"lifetime_seconds" in the ``mechanism_opts`` dict.  Likewise, the second option, ``'CHEMISTRY/initial_cond/0/concentration'``
has keys for the first, second, and fourth parts of the path , but the ``initial_cond`` value is a list, so the third index is 
the numeric list index ``0``.

The second argument to the ``add_ens_var_by_string`` function is the values that the option should have in each of the ensemble members.
This is where the ``ensemble_mode`` argument for the ``EnsembleRunner`` comes it. It can have two values:

- ``'iterations'`` means that each option modified must have the same number of values, :math:`n`. The ensemble will have :math:`n`
  members, and for each member :math:`i` (where :math:`0 \leq i < n`), the modified options will have the value at index :math:`i`.
- ``'combinations'`` means that each option modified can have any number of values, and the ensemble will consist of all possible 
  combinations of those values.

In our case, we chose ``'iterations'``, so our ensemble will have three members:

======  ========  ===================
Member  Lifetime  Init. concentration
======  ========  ===================
1       3600      2
2       7200      20
3       10,800    200
======  ========  ===================

If instead we had chosen ``'combinations'``, our ensemble would have nine members:

======  ========  ===================
Member  Lifetime  Init. concentration
======  ========  ===================
1       3600      2
2       3600      20
3       3600      200
4       7200      2
5       7200      20
6       7200      200
7       10,800    2
8       10,800    20
9       10,800    200
======  ========  ===================


For each of these members, the ensemble runner will create a new directory in :file:`test_output` to write to.
By default, these directories will have the name ``pecans_ens_member_INDEX``, with INDEX being the ensemble member
number (starting from 0). You can change this - see the ``member_naming_fxn`` argument of :class:`~pecans.ensembles.api.EnsembleRunner`.
Each of those directories will have all the output files from its respective ensemble member's run.

There's one last step, we have to execute the ensemble by calling its ``run`` method. Put all together, this example is::

    from pecans.ensembles.api import EnsembleRunner

    runner = EnsembleRunner(
        base_config_file = 'examples/one_box_ideal/one_box_ideal.toml',
        ensemble_mode='iterations',
        root_output_dir='test_output'
    )

    runner.add_ens_var_by_string('CHEMISTRY/mechanism_opts/lifetime_seconds', [3600, 7200, 10800])
    runner.add_ens_var_by_string('CHEMISTRY/initial_cond/0/concentration', [2, 20, 200])

    runner.run()

Ensemble API functions
----------------------

.. include:: includes/top_includes.rst

.. automodule:: pecans.ensembles.api
   :members: