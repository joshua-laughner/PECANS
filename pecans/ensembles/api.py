from collections import OrderedDict
import copy
import itertools
import os
import warnings

from ..utilities.config import load_config_file
from ..main import Domain

from typing import Optional, Sequence, Any


def _make_member_name(member_index, **config_opts):
    return 'pecans_ens_member_{:06}'.format(member_index)


class EnsembleError(Exception):
    """Error type used for problems in setting up an ensemble run
    """
    pass


class EnsembleRunner:
    """
    Manager class to help run an ensemble of PECANS instances.

    :param base_config_file: the path to the configuration file to use as a starting point for the ensemble.

        .. note::
           The ``output_path`` argument in the configuration is ignored; the output path is set by the 
           ``root_output_dir``,  ``member_naming_fxn``, and ``save_in_individual_dirs`` options of this
           class.

    :param ensemble_variables: the configuration variables to modify for the different ensemble members. This must be
        a dictionary where the keys are the options given as strings and the values are lists of the values that each
        ensemble member will have (see the :ref:`tutorial <ensemble_tutorial>` for details). The keys will have the 
        form ``"key1/key2/key3"`` and map to the configuration as ``config[key1][key2][key3]``. For list configuration 
        elements, give the index as part of the key (e.g. ``"CHEMISTRY/initial_cond/0/concentration"``) and it will 
        be automatically converted to an integer if needed.

        Alternatively, you can construct the ensemble runner without this argument and add values later using the
        :meth:`~pecans.ensembles.api.EnsembleRunner.add_ens_var_by_string` method.

    :param ensemble_mode: this determines how the ensemble_variables are varied. Possible options are:

        * ``'iterations'`` - each ensemble variable has its values specified in an iterable (i.e. list or numpy array); the
          nth ensemble member will use the nth value for each variable.
        * ``'combinations'`` - each ensemble variable has its possible values specified in an iterable, but unlike
          ``'iterations'``, all possible combinations are tested. 
        
        If running in ``'iterations'`` mode, then the length of the value lists in ``ensemble_variables`` must all be
        equal.

        See the :ref:`tutorial <ensemble_tutorial>` for detailed examples.

    :param root_output_dir: the root directory to place the ensemble output in. Default is the current directory.

    :param member_naming_fxn: a function that accepts the ensemble member index as an integer and the member's options 
        as keyword arguments and returns a string that should be a unique name for that ensemble member. The default will 
        return ``"pecans_ens_member_N"`` where N is the ensemble member index. This name will be used for the member's output
        directory name (if ``save_in_individual_dirs`` is ``True``) and, with ".nc" appended, the output file name if 
        ``save_final_output_only`` is True).

        You can use this to set the output names to something that incorporates the varied ensemble variables into the
        file/directory names. For example, if dx and dy are being varied, you could do::

            def custom_name(member_index, **config_opts):
                dx = config_opts['DOMAIN/dx'] / 1000  # convert to kilometers
                dy = config_opts['DOMAIN/dy'] / 1000

                return 'pecans_ens_dx-{}km_dy-{}km'.format(dx, dy)

            ensemble = EnsembleRunner( ... , member_naming_fxn=custom_name)

        This would put the dx and dy values (converted to kilometers) into the file or directory names.

    :param save_in_individual_dirs: optional, default is ``True``. This creates separate directories for the output of each
        ensemble member, named using the member_naming_fxn. Each member's output is saved in the corresponding directory.

    :param save_final_output_only: optional, default is ``False``, meaning that each ensemble member will save output at the
        output frequency defined in their configuration. If ``True``, only the final state of the model will be saved, and
        it will be named by the member_naming_fxn plus the '.nc' extension. If this is ``False``, save_in_individual_dirs
        will automatically be set to ``True``; save_in_individual_dirs is not already ``True``, a warning is issued.
    """
    @property
    def base_config(self):
        return self._base_config

    @property
    def ensemble_mode(self):
        return self._mode

    def __init__(self, base_config_file: str, ensemble_variables: Optional[dict] = None, ensemble_mode: str = 'iterations',
                 save_in_individual_dirs: bool =True, save_final_output_only: bool =False, member_naming_fxn=_make_member_name,
                 root_output_dir: Optional[str] = '.'):

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
            config['OUTPUT']['output_frequency'] = 0
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

    def add_ens_var_by_string(self, ensemble_variable: str, values: Sequence[Any]):
        """
        Add an option that should be varied among the different ensemble members.

        :param ensemble_variable: The option to vary as a string, must follow the same format as the keys for the
            ``ensemble_variables`` argument of the class constructor.

        :param values: values describing how that option should be varied among the different ensemble members. The
            required form varies depending on the ensemble_mode option set during initialization. See the class
            documentation for specifics.

        :return: none
        """

        self._validate_ensemble_variable_string(self.base_config, ensemble_variable)
        self._ensemble_variables[ensemble_variable] = values

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
            member_config['OUTPUT']['output_frequency'] = 0

        member_domain = Domain(member_config, output_dir=output_dir)
        member_domain.execute()

        if self._final_output_only:
            member_domain.write_output(os.path.join(self._root_output_dir, member_out_name + '.nc'))

    @classmethod
    def _modify_member_config(cls, config, config_opts):
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
            # First we need to recursively descend through the configuration to find the element to modify
            # Since the configuration is a TOML file, all the intermediate levels should be dicts or lists,
            # so we can interpret our index based on the element type.
            curr_config_element = config
            curr_config_path = []

            opt_parts = opt.split('/')
            for ipart, part in enumerate(opt_parts):
                if isinstance(curr_config_element, dict):
                    idx = part
                elif isinstance(curr_config_element, list):
                    try:
                        idx = int(part)
                    except ValueError:
                        path = '/'.join(curr_config_path)
                        raise EnsembleError(f'Configuration element {path} is a list, so the next part of the path must be an integer, but got {part}')
                else:
                    path = '/'.join(curr_config_path)
                    raise EnsembleError(f'{path} is not a dictionary, list, or tuple in the configuration (it is a {type(curr_config_element)}), so it cannot be indexed by {part}')

                # Don't actually do the last indexing because we need to modify the value in its parent element
                if ipart < len(opt_parts) - 1:
                    try:
                        curr_config_element = curr_config_element[idx]
                    except (IndexError, KeyError) as e:
                        path = '/'.join(curr_config_path)
                        raise EnsembleError(f'{idx} is not a valid index for {path}: {e}')
                elif isinstance(curr_config_element, dict) and idx not in curr_config_element:
                    # Test that the key exists for dictionaries so that we're not accidentally creating a new key
                    path = '/'.join(curr_config_path)
                    raise EnsembleError(f'Cannot set value of {part} in {path}, "{part}" is not a key in that section')
                else:
                    # Last part of the path, so we set the value
                    try:
                        curr_config_element[idx] = val
                    except IndexError:
                        # Easier to handle list index checks here
                        raise EnsembleError(f'Cannot set value of index {idx} in {path}, {idx} is not a valid index of that list')
                curr_config_path.append(part)
                        
        

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
        
        # To validate we'll just test if we can modify the configuration - that method has the error checking built in
        config = copy.deepcopy(config)
        cls._modify_member_config(config, {ensemble_var: None})
        

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

        if section not in config:
            error_helper('No section "{section}" in the configuration'.format(section=section))
        elif option not in config[section]:
            error_helper('No option "{option}" in section ("{section}")in the configuration'.format(option=option, section=section))
        elif suboption is not None:
            # If there's suboption, then we need to check that the specified option is a dictionary
            # (or list eventually) and that it has that key
            config_opt = config[section][option]
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
