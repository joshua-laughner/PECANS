Using PECANS
============

Installation
------------

Downloading PECANS
~~~~~~~~~~~~~~~~~~

We strongly recommend you clone the PECANS repository to get the latest version. To do so, first make sure
you have Git installed on your computer by running ``git --version`` in a terminal. As long as it returns
something like ``git version 2.39.3``, you have it. If you do not, search online for instructions to install
Git for your computer type (Windows, Mac, or Linux).

Then, in terminal go to the directory where you want to download PECANS and run::

    git clone https://github.com/joshua-laughner/PECANS.git

This will create the PECANS directory in your current directory.


Getting the correct Python version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PECANS is currently built for Python 3.11. Other Python versions may work, but will not be supported.
If you use [conda](https://docs.conda.io/en/latest/) as a package manager (often because you installed
[Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)),
you can install an isolated Python 3.11 with the command::

    conda create --name py311 python=3.11

The important part is the last argument, ``python=3.11``. This tells conda to install Python 3.11 in your
new environment. You can use whatever name you wish for the environment; we use ``py311`` here because this
environment will only have Python 3.11.

There are many other tools you could use to install Python 3.11. Whatever approach you take, make sure that
this Python 3.11 is your active Python version for the next step. If you followed the conda command above,
you can make this so with the command::

    conda activate py311

You can confirm this worked by the output of ``python --version``; if it is some version of 3.11, you are set.
Note that conda environments are only activated in your current shell - so if you switch to another terminal 
tab or window (or close and reopen terminal), this will *not* be you active python version in that shell.

.. _installing_dependencies:

Installing dependencies
~~~~~~~~~~~~~~~~~~~~~~~

For now, PECANS uses the old-style pip requirements file to manage dependencies. We recommend installing these
dependencies in a virtual environment to avoid conflicts with any existing Python tools you have.

First, make sure that Python 3.11 is your active Python version with the ``python --version``. If the output is similar
to ``Python 3.11.5``, then you are ready to proceed. (The last number, "5" here, may be different for you.)

Second, navigate to the top directory of the PECANS repository you downloaded from GitHub. That will be the directory
containing the ``run_pecans.py`` and ``requirements.txt`` files, among others.

Third, create a virtual environment with the command ``python -m venv .venv``. This will create a new hidden directory,
``.venv``, in the PECANS repo that we will install our dependencies into.

Fourth, activate this virtual environment. We need to know what shell we are using, so first run ``echo $SHELL`` and note
what it prints out. Now, from the same directory as the previous step (the PECANS repo root), run *one* of the following
commands:

- If ``echo $SHELL`` prints out a path ending in ``bash`` or ``zsh``, run ``source .venv/bin/activate``
- If ``echo $SHELL`` prints out a path ending in ``csh`` or ``tcsh``, run ``source .venv/bin/activate.csh``
- If ``echo $SHELL`` prints out something else, I'm assuming you installed your own shell and don't need me to tell you what to do.

Fifth, ensure you have the latest version of ``pip`` by running ``pip install --upgrade pip``. This will have the ``pip`` installed in
our virtual environment upgrade itself if needed.

Sixth, install PECANS' dependencies by running ``pip install -r requirements.txt``. This will install all the packages listed in the
``requirements.txt`` file into the virtual environment.

.. _install_pecans_package:

Install PECANS as a package (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While PECANS is primarily meant to be run from the command line, you may also want to be able to import parts of it into other Python
programs or Jupyter notebooks. To allow that, activate the PECANS virtual environment (step four in :ref:`installing_dependencies`),
then run the command ``python setup.py develop`` from the root of the PECANS repo. This will install PECANS as an "editable" package into the
same virtual environment as its dependencies. Making it an "editable" package means that any changes to the code will take effect the
next time you import PECANS into your projects, without needing to rerun this ``pip`` command.

.. note::
   You will probably get a warning about the setup.py install being deprecated. We're aware of this, but using the recommended approach
   (``pip install -e .``) does not currently work with Cython, which is needed for the compiled mechanisms.

.. note::
   To import PECANS into Python scripts or notebooks, the script/notebook must be run with the virtual environment. For scripts, that
   is just a matter of activating the virtual environment before running your script with ``python YOUR_SCRIPT``.  For Jupyter notebooks,
   that means adding a kernel that knows to run in the virtual environment. To do that, I normally follow 
   `this approach <https://stackoverflow.com/a/53546675>`_. Since we're working with a vanilla virtual environment rather than a conda environment, 
   the commands are slightly different than given in that link. With the virtual environment active, run ``pip install ipykernel``
   then the ``ipython kernel install --user --name=<any_name_for_kernel>`` command. The next time you start Jupyter, you should see this
   kernel available under the name you gave to the ``--name`` option.


Running PECANS
--------------

For all the steps here, make sure you have the virtual environment we created :ref:`above <installing_dependencies>` active. See the Fourth
step in :ref:`installing_dependencies` if you need a reminded how to activate it.

First, you must create a configuration file for your run. There are several examples in the :file:`examples` directory of the repo.
Your configuration file must be ``pecans_config.toml`` in the root of the PECANS repo. (That can be a symbolic link to another configuration
file - try linking the :file:`examples/one_box_ideal/one_box_ideal.toml` file for your first run.) The configuration file is covered in
detail under :ref:`config_user_guide`, so read that section and check out the examples to learn how to create these files.

Second, if you want to use a realistic chemical mechanism, rather than one of the :ref:`idealized <ideal_chem_mech>` ones, it will need
built. There are several demo mechanisms included with PECANS. To build one, run the :file:`build_pecans.py` script (``python build_pecans.py``,
assuming you're in the PECANS repo root directory). This will present a list of the available mechanisms. For your first run, choose the
"nox" mechanism. The build script should run rather quickly, and when it finished, you should see a file named like :file:`chemderiv.*.so`
in the :file:`pecans` directory. (If you only see :file:`chemderiv.pyx`, then the mechanism did not build fully - see the note for advice in this case.)

.. note::
   Although :file:`build_pecans.py` is set up to automatically compile the mechanism to the :file:`.so` library, sometimes Cython does not
   cooperate. If running :file:`build_pecans.py` doesn't produces the :file:`chemderiv.pyx` file but not the :file:`chemderiv.*.so` file, try
   installing PECANS as a package into the virtual environment as described in :ref:`install_pecans_package` and rerunning :file:`build_pecans.py`. 
   Having PECANS installed as a package seems to make Cython behave a little better.

Finally, all that you need to do to run PECANS is call ``python run_pecans.py``. If you linked the ``one_box_ideal`` example configuration, you will 
see some output files produced in the :file:`test_output` directory. You can use the :file:`plot_pecans.py` script for some quick plots of them, or
read the files yourself and check out the results.