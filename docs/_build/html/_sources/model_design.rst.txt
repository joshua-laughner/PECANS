Model Design
============

At its core, the PECANS model is built around the :class:`~pecans.core.Domain` object, which contains the model
concentrations within it and is linked to various solvers that handle solving the necessary differential equations for
a single physical process. For example, right now, there are solvers for chemistry, transport and emissions. By
`Strang operator splitting <https://en.wikipedia.org/wiki/Strang_splitting>`_, each of these solvers can be applied in
turn during a time step, and so in the :func:`~pecans.core.Domain.step` method, each solver is applied to the array of
concentrations separately.

PECANS allows for the possibility of different solvers, for example, some methods are more accurate but take more time
to compute, so for larger models, you may find it necessary to switch to a cheaper, less accurate method. To allow for
this flexibility, this means that the solvers need to be set up first.


Design philosophy
-----------------

Anyone who wants to contribute to this project is welcome to! But please read this section and take it seriously. It is
important that this model adhere to certain practices to help make the code easily readable and understood.

#. **Constant globals only**: Data and variables should be explicitly passed into any functions that need them. Relying
   on global/module variables to carry mutable data is confusing and difficult to follow. Constant values used
   throughout a module can and should be defined as module variables, but these must not be changed.

#. **No use of from, import \***: This makes it unclear where variables or functions came from, since there is no
   explicit indication in this import statement.

Model configuration
-------------------

The specifics of the model are determined by the `pecans_config.cfg` file in the same directory as `run_pecans.py`. At
the beginning of a model run, this is ingested and represented by a :class:`~pecans.utilities.config.BetterConfig`
instance. This is derived from :class:`configparser.RawConfigParser`, and so behaves similarly in that options are
organized into sections, but differently in that option values are automatically parsed into Python literals if possible,
and only kept as a string if its type cannot be inferred. This instance will need to be passed to any functions that
need it.


Solver organization
-------------------

Each solver is placed in a sub-package of the main PECANS package, e.g. chemistry is handled by the chemistry package,
transport by the transport package, etc. Within each package should be a <package>_setup module, (e.g. chemistry_setup,
transport_setup) that has the setup_<package> method (e.g. setup_chemistry(), setup_transport(), etc.). These setup
methods should require only one argument, the configuration object.

This setup function needs to return the driver function that will actually solve the differential equations. If these
solvers rely on information other than the chemical species' concentrations that may change every timestep, that data
should be read in within the driver function (e.g. the transport driver should call the function to read the current
times' wind fields within it, rather than rely on this being passed in from the domain). Different drivers may be
defined that solve the differential equations differently, but all driver functions must accept the same inputs and
return the same outputs. (Use of the ``**kwargs`` construct will be useful to consume extra keyword arguments not used.)
When called from the domain object, these driver functions should have all inputs specified by keyword, not position.
That way if a new driver is added that requires additional inputs, it is easy to add those to the call in the Domain
object.