from collections import OrderedDict
import copy
import itertools
import os
import warnings

from ..utilities.config import load_config_file
from ..main import Domain

import pdb

# TODO:
#   Document 'iterations' and 'combinations' options and the ensemble runner in general


def _make_member_name(member_index, **config_opts):
    return 'pecans_ens_member_{:06}'.format(member_index)


class EnsembleRunner:
    """
    Manager class to help run an ensemble of PECANS instances.

    :param base_config_file: the path to the configuration file to use as a starting point for the ensemble.
    :type base_config_file: str

    :param ensemble_variables: the configuration variables to modify for the different ensemble members. This must be
        a dictionary where the keys are the options given as 'SECTION/OPTION' (e.g. 'DOMAIN/dx' to set the dx value)
        or 'SECTION/OPTION/SUBOPTION', if the option contains a suboption (e.g. 'CHEMISTRY/mechanism_opts/lifetime_seconds'
        to set the lifetime for the ideal_first_order mechanism), and the values describe how to vary that option. The
        exact format that the values take on depends on the ensemble mode (see next parameter).

        ensemble_variables does not need to be specified during initialization of this class; additional variables can
        be specified using the add_ens_var_by_string() or add_ens_var_by_option() methods, if that is more convenient.
    :type ensemble_variables: dict or None

    :param ensemble_mode: this determines how the ensemble_variables are varied. Possible options are:

        * 'iterations' - each ensemble variable has its values specified in an iterable (i.e. list or numpy array); the
          nth ensemble member will use the nth value for each variable.
        * 'combinations' - each ensemble variable has its possible values specified in an iterable, but unlike
          'iterations', all possible combinations are tested. E.g. if your ensemble variables are
          {'DOMAIN/dx': [1000, 4000, 9000], 'DOMAIN/dy': [2000, 5000, 8000]} then a total of 9 ensemble members will be
          run, one with dx = 1000 & dy = 2000, one with dx = 1000 & dy = 5000, one with dx = 4000 & dy = 2000, etc.

    :type ensemble_mode: str

    :param root_output_dir: the root directory to place the ensemble output in. Default is the current directory.
    :type root_output_dir: str

    :param member_naming_fxn: a function that accepts two arguments (the ensemble member index as an integer and the
        member's options as keyword arguments) and returns a string that should be a unique name for that ensemble
        member. The default will return "pecans_ens_member_N" where N is the ensemble member index. This name will be
        used for the member's output directory name (if save_in_individual_dirs is True) and, with ".nc" appended, the
        output file name if save_final_output_only is True).

        You can use this to set the output names to something that incorporates the varied ensemble variables into the
        file/directory names. For example, if dx and dy are being varied, you could do::

            def custom_name(member_index, **config_opts):
                dx = config_opts['DOMAIN/dx'] / 1000  # convert to kilometers
                dy = config_opts['DOMAIN/dy'] / 1000

                return 'pecans_ens_dx-{}km_dy-{}km'.format(dx, dy)

            ensemble = EnsembleRunner( ... , member_naming_fxn=custom_name)

        This would put the dx and dy values (converted to kilometers) into the file or directory names.
    :type member_naming_fxn: callable

    :param save_in_individual_dirs: optional, default is True. This creates separate directories for the output of each
        ensemble member, named using the member_naming_fxn. Each member's output is saved in the corresponding
        directory.
    :type save_in_individual_dirs: bool

    :param save_final_output_only: optional, default is False, meaning that each ensemble member will save output at the
        output frequency defined in their configuration. If True, only the final state of the model will be saved, and
        it will be named by the member_naming_fxn plus the '.nc' extension. If this is False, save_in_individual_dirs
        will automatically be set to True; save_in_individual_dirs is not already True, a warning is issued.
    :type save_final_output_only: bool
    """
    @property
    def base_config(self):
        return self._base_config

    @property
    def ensemble_mode(self):
        return self._mode

    def __init__(self, base_config_file, ensemble_variables=None, ensemble_mode='iterations',
                 save_in_individual_dirs=True, save_final_output_only=False, member_naming_fxn=_make_member_name,
                 root_output_dir='.'):

        config = load_config_file(base_config_file)
        self._base_config = config

        # Store ensemble_variables to an ordered dictionary since although the input order probably doesn't matter,
        # having the variables in a consistent order can help with some of the ensemble modes.
        self._ensemble_variables = OrderedDict()
        if ensemble_variables is not None:
            self._validate_init_ensemble_variables(config, ensemble_variables)

            for key, val in ensemble_variables.items():
                self.add_ens_var_by_string(key, val)

        # If save_final_output_only is False and save_in_individual_dirs is False, there's going to be issues with
        # overwriting files. We'll force save_in_individual_dirs to be False if save_final_output_only is not, but
        # issue a warning that that's what we're doing.
        if not save_final_output_only:
            if not save_in_individual_dirs:
                save_in_individual_dirs = True
                warnings.warn('With save_final_output_only = False, save_in_individual_dirs must be True to avoid '
                              'different ensemble members overwriting each other\'s output. Changing '
                              'save_in_individual_dirs to True.')

        self._root_output_dir = root_output_dir
        self._member_naming_fxn = member_naming_fxn
        self._save_indiv_dirs = save_in_individual_dirs
        self._final_output_only = save_final_output_only

        if save_final_output_only:
            config.set('OUTPUT', 'output_frequency', 0)
        elif not save_in_individual_dirs:
            raise RuntimeError('If not set to save the final output only (save_final_output_only=False), '
                               'the ensemble must be set to save model output in individual directories '
                               '(save_in_individual_dirs=True)')

        if ensemble_mode == 'iterations':
            self._run_fxn = self._iterations_driver
        elif ensemble_mode == 'combinations':
            self._run_fxn = self._combinations_driver
        else:
            raise NotImplementedError('No run function implemented for ensemble_mode={}'.format(ensemble_mode))

    def add_ens_var_by_string(self, ensemble_variable, values):
        """
        Add an option that should be varied among the different ensemble members.

        :param ensemble_variable: The option to vary as a string, 'SECTION/OPTION' (e.g. 'DOMAIN/dx') or
            'SECTION/OPTION/SUBOPTION' (e.g. 'CHEMISTRY/mechanism_opts/lifetime_seconds').
        :type ensemble_variable: str

        :param values: values describing how that option should be varied among the different ensemble members. The
            required form varies depending on the ensemble_mode option set during initialization. See the class
            documentation for specifics.

        :return: none
        """

        self._validate_ensemble_variable_string(self.base_config, ensemble_variable)
        self._ensemble_variables[ensemble_variable] = values

    def add_ens_var_by_option(self, section, option, values_or_subopt, subopt_values=None):
        """
        An alternate method to add an option that should be varied among the different ensemble members.

        This method will be called with three or four arguments, depending on if there is a suboption that needs to be
        specified. If not::

            add_ens_var_by_option( section, option, values)

        If so::

            add_ens_var_by_option( section, option, suboption, values )

        This is just a convenience form of add_ens_var_by_string that allows you to specify the section, option, and
        suboption separately instead of combining them into a string.

        :return: none
        """
        # Allow for the two or three argument forms (section, option, values) and (section, option, suboption, values)
        if subopt_values is None:
            values = values_or_subopt
            suboption = None
        else:
            values = subopt_values
            suboption = values_or_subopt

        opt_string = '{}/{}'.format(section, option)
        if suboption is not None:
            opt_string += '/{}'.format(suboption)
        self.add_ens_var_by_string(opt_string, values)

    def run(self):
        """
        Carry out all the ensemble simulations.

        :return: none
        """
        self._run_fxn()

    def run_one_member(self, member_index, config_opts):
        """
        Run a single member of the ensemble.

        :param member_index: a unique index identifying which member of the ensemble this is.
        :type member_index: int

        :param config_opts: dictionary where the keyword specifies the option to change from the base
            configuration. For example::

                run_one_member(0, {'DOMAIN/dx':1000, 'DOMAIN/dy':2000})

            will run a member of the ensemble (#0) with the dx value set to 1000 m and dy to 2000 m.

        :return: none

        This method is usually called internally by the ``run`` method, which automatically iterates over all ensemble
        members, but in some cases you may want more control over how the different members are run, while still taking
        advantage of this method's built-in capability to only modify a few of the configuration options and redirect
        model output to a member-specific file or directory.
        """
        member_config = copy.deepcopy(self.base_config)
        self._modify_member_config(member_config, config_opts)
        member_out_name = self._member_naming_fxn(member_index, **config_opts)
        if self._save_indiv_dirs:
            output_dir = os.path.join(self._root_output_dir, member_out_name)
            os.makedirs(output_dir, exist_ok=True)
        else:
            output_dir = self._root_output_dir
            if not os.path.isdir(output_dir):
                raise IOError('Specified output directory ({}) does not exist'.format(output_dir))

        if self._final_output_only:
            # Setting the output frequency to 0 disables automatic writing of output by the domain. This way we can
            # manually save the final output with a more meaningful name
            member_config.set('OUTPUT', 'output_frequency', 0)

        member_domain = Domain(member_config, output_dir=output_dir)
        member_domain.execute()

        if self._final_output_only:
            member_domain.write_output(os.path.join(self._root_output_dir, member_out_name + '.nc'))

    def _modify_member_config(self, config, config_opts):
        """
        Internal helper method that handles modifying certain options in a given configuration.

        :param config: the configuration instance to modify.
        :type config: :class:`~pecans.utilities.config.BetterConfig`

        :param config_opts: dictionary specifying the SECTION/OPTION or SECTION/OPTION/SUBOPTION to modify and
            the value to set it to, e.g. ``_modify_member_config(config, {'DOMAIN/dx':1000, 'DOMAIN/dy':2000}``.
        :type config_opts: dict

        :return: none
        """
        for opt, val in config_opts.items():
            opt_parts = opt.split('/')
            if len(opt_parts) == 2:
                section, option = opt_parts
                config.set(section, option, val)
            elif len(opt_parts) == 3:
                section, option, suboption = opt_parts
                config.get(section, option)[suboption] = val
            else:
                raise NotImplementedError('Cannot handle config option "{}": expected 2 or 3 parts separated by slashes'.format(opt))

    def _run_ensemble(self, config_opts):
        """
        Internal function that handles running the entire ensemble.

        :param config_opts: dictionary specifying the SECTION/OPTION or SECTION/OPTION/SUBOPTION to modify and
            the values each one should be set to for each member, as an iterable, e.g.
            ``_run_ensemble({'DOMAIN/dx': [1000,4000,9000], 'DOMAIN/dy': [2000,5000,8000]}`` would run a three member
            ensemble where the first member has dx = 1000 & dy = 2000, the second dx = 4000 & dy = 5000, and the third
            dx = 9000 & dy = 8000. The length of each value in the dictionary must be the same.
        :type config_opts: dict

        :return: none
        """
        base_key = [k for k in config_opts.keys()][0]
        opt_length = len(config_opts[base_key])
        for key, val in config_opts.items():
            if len(val) != opt_length:
                raise ValueError('The number of values given for "{}" is inconsistent with those given for "{}"'.format(
                    key, base_key
                ))

        # Later this can be modified to run in parallel
        for i_ens in range(opt_length):
            member_opts = {k: v[i_ens] for k, v in config_opts.items()}
            self.run_one_member(i_ens, member_opts)

    def _iterations_driver(self):
        """
        Driver function for ensemble_mode = 'iterations'

        :return: none
        """
        self._run_ensemble(self._ensemble_variables)

    def _combinations_driver(self):
        """
        Driver function for ensemble_mode = 'combinations'

        :return: none
        """
        # We need to create a new version of the ensemble opts dictionary that has all possible combinations of each of
        # the input variables. We take advantage of the ensemble_variables being converted to an OrderedDict while doing
        # this, so that we are sure that the keys and values will always present in the same order.
        #
        # The itertools.product() function will return lists that are each possible combination of elements from all the
        # input lists. Since the values and keys are ordered, this means than combo[0] corresponds to key[0] and so on.
        iter_ens_vars = {k: [] for k in self._ensemble_variables.keys()}
        for combo in itertools.product(*[v for v in self._ensemble_variables.values()]):
            for i_key, key in enumerate(self._ensemble_variables.keys()):
                iter_ens_vars[key].append(combo[i_key])

        self._run_ensemble(iter_ens_vars)

    @classmethod
    def _validate_init_ensemble_variables(cls, config, ensemble_vars):
        """
        Internal helper function that checks that the ensemble variables given during initialization are valid

        :param config: the base configuration object read in, variables are checked to ensure they exist in the config.
        :type config: :class:`~pecans.utilities.config.BetterConfig`

        :param ensemble_vars: the dictionary of ensemble variables given to the __init__ method.
        :type ensemble_vars: dict

        :return: none, raises ValueError or TypeError if ensemble_vars is incorrectly formatted.
        """
        if not isinstance(ensemble_vars, dict):
            raise TypeError('ensemble_variables must be a dictionary')

        # Next, is every key a slash separated string, where the elements separated by slashes are a valid section,
        # option, and (if necessary) suboption
        for key in ensemble_vars.keys():
            cls._validate_ensemble_variable_string(config, key)

    @classmethod
    def _validate_ensemble_variable_string(cls, config, ensemble_var):
        """
        Internal helper function that validates a single ensemble option string

        :param config: the base configuration object read in, variables are checked to ensure they exist in the config.
        :type config: :class:`~pecans.utilities.config.BetterConfig`

        :param ensemble_var: the string form of the option name, e.g. 'DOMAIN/dx' or
            'CHEMISTRY/mechanism_opts/lifetime_seconds'.
        :type ensemble_var: str

        :return: none

        Internally, the string is parsed and the section, option, and suboption passed to _validate_ensemble_variable
        to check that they exist in the configuration.
        """
        if not isinstance(ensemble_var, str):
            raise TypeError('_validate_ensemble_variable_string expects the ensemble variable as a string')
        elif ensemble_var.count('/') < 1:
            raise ValueError('An ensemble variable specified as a string must contain at least the section and option'
                             ' to modified separated by a slash, e.g. DOMAIN/nx.')
        else:
            key_elements = ensemble_var.split('/')
            if len(key_elements) == 2:
                section, option = key_elements[0:2]
                cls._validate_ensemble_variable(config, section, option, original_key=ensemble_var)
            elif len(key_elements) == 3:
                section, option, suboption = key_elements
                cls._validate_ensemble_variable(config, section, option, suboption, original_key=ensemble_var)
            else:
                raise ValueError('An ensemble variable specified as a string may have at most a section, option, '
                                 'and suboption (this error is raised if > 2 slashes exist in the string')

    @staticmethod
    def _validate_ensemble_variable(config, section, option, suboption=None, original_key=None):
        """
        Internal helper function that checks that the requested section, option, suboption exist in the configuration.

        :param config: the base configuration object read in, variables are checked to ensure they exist in the config.
        :type config: :class:`~pecans.utilities.config.BetterConfig`

        :param section: the section name in the configuration
        :type section: str

        :param option: the option name in the configuration
        :type option: str

        :param suboption: the suboption name, must be omitted or set to None if the option does not have suboptions.
        :type suboption: str or None

        :param original_key: the original full string of the option name, e.g. 'DOMAIN/dx'. Optional, just used to
            provide better error messages if an improperly formatted option was parsed from a string.

        :return: none, raises ValueError if any error with the option specification.
        """
        def error_helper(msg):
            if original_key is not None:
                # If the caller provided the original key (i.e. the string parsed, like "CHEMISTRY/do_chemistry"),
                # then we can add that to the error message to be more helpful.
                #
                # Remove any terminating period from the error message before adding the additional information.
                msg = msg.rstrip('.') + ' (from ensemble variable "{key}").'.format(key=original_key)
            raise ValueError(msg)

        if not config.has_section(section):
            error_helper('No section "{section}" in the configuration'.format(section=section))
        elif not config.has_option(section, option):
            error_helper('No option "{option}" in section ("{section}")in the configuration'.format(option=option, section=section))
        elif suboption is not None:
            # If there's suboption, then we need to check that the specified option is a dictionary
            # (or list eventually) and that it has that key
            config_opt = config.get(section, option)
            if not isinstance(config_opt, dict):
                error_helper(
                    'A suboption ({subopt}) was given, but the option "{option}" does not have suboptions'.format(
                        subopt=suboption, option=option
                    ))
            elif suboption not in config_opt:
                error_helper('The suboption "{subopt}" is not a suboption specified in the {section}/{option} '
                                 'option.'.format(subopt=suboption, section=section, option=option))

    @staticmethod
    def _validate_ensemble_opts(ensemble_mode, ensemble_opts):
        """
        Internal helper function that ensures that the required ensemble options are present for a given mode

        :param ensemble_mode: the ensemble mode to check options for
        :type ensemble_mode: str

        :param ensemble_opts: dictionary of keyword options passed to the __init__ method
        :type ensemble_opts: dict

        :return: none, raises ValueError if option missing or otherwise incorrect
        """
        # This dictionary of dictionary defines the required options for each ensemble mode; within each child
        # dictionary, the keys give the keyword and the set value gives the allowed values that option may take on. An
        # empty set indicates that any value is permitted; functions may also be given

        ensemble_opt_requirements = {
            'iterations': dict(),
            'combinations': dict(),
        }

        if ensemble_mode not in ensemble_opt_requirements.keys():
            raise ValueError('"{}" is not a recognized ensemble mode'.format(ensemble_mode))

        requirements = ensemble_opt_requirements[ensemble_mode]
        for opt, allowed_vals in requirements.items():
            if opt not in ensemble_opts:
                raise ValueError('With ensemble_mode = {mode}, "{opt}" is a required ensemble option'.format(
                    mode=ensemble_mode, opt=opt
                ))

            err_msg = '{val} is not a legal value for the "{opt}" option.'.format(val=ensemble_opts[opt], opt=opt)
            if isinstance(allowed_vals, set):
                if ensemble_opts[opt] not in allowed_vals and len(allowed_vals) > 0:
                    raise ValueError(err_msg)
            else:
                if not allowed_vals[ensemble_opts[opt]]:
                    raise ValueError(err_msg)
