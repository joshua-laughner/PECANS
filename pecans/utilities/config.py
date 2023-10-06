from ast import literal_eval
import configparser
import copy
import os
import re
import tomllib

import pdb

# Adapted from https://stackoverflow.com/a/26634150, this regex can be used to split on commas that are outside of
# brackets, i.e. [ ], only. It looks to find a comma that is NOT FOLLOWED by any number of non-bracket characters then a
# closing bracket.
# split_on_comma_re = re.compile(r',(?![^\[\]]*\])')
# This will be useful in the future, for now we're just going to split on commas.
split_on_comma_re = re.compile(r',')


class ConfigurationError(Exception):
    """
    An error that represents something wrong in the PECANS configuration.

    This should be raised specifically if there is something wrong with what the user has given in the config file.
    """
    pass


def load_config_file(config_file):
    """
    Reads a configuration file and returns a BetterConfig instance representing it

    This function automatically instantiates a :class:`~pecans.utilities.config.BetterConfig` instance, sets it to parse
    values into Python literals (including tuples or dicts)
    :param config_file: the file to read.
    :type config_file: str or list
    :return: the :class:`~pecans.utilities.config.BetterConfig` object
    """
    with open(config_file, 'rb') as f:
        return tomllib.load(f)


def get_domain_size_from_config(config, all_dims=False):
    """
    Helper function that reads the domain size (as a tuple) from a configuration object.

    :param config: the configuration object to read the domain size from.
    :type config: :class:`~pecans.utilities.config.BetterConfig`

    :param all_dims: whether to include just non-zero dimensions (False, default) or all dimensions (True)
    :type all_dims: bool

    :return: the shape of the domain, as a tuple with 1, 2, or 3 elements depending if the domain is 1-, 2-, or 3- D.
    :rtype: tuple of int
    """
    ny = config['DOMAIN']['ny']
    nx = config['DOMAIN']['nx']
    nz = config['DOMAIN']['nz']

    # Check that if nz > 0, so are nx and ny; if ny > 0, so is nx, because, e.g. nx, nz > 0 but ny == 0 is an
    # invalid setup (2D must be along x and y).
    if nz > 0:
        if nx <= 0 or ny <= 0:
            raise ConfigurationError('If nz > 0, nx and ny must also both be > 0')
        else:
            shape = (nx, ny, nz)
    elif ny > 0:
        if nx <= 0:
            raise ConfigurationError('If ny > 0, nx must also be > 0')
        else:
            shape = (nx, ny)
    elif nx < 1:
        raise ConfigurationError('nx must be >= 1')
    else:
        shape = (nx,)

    # Include the above if statements to do the error checking
    if all_dims:
        shape = (nx, ny, nz)

    return shape


def list_missing_opts(required_opts, config, section, raise_error=False):
    section_dict = config[section]
    missing_opts = _opts_missing_from_dict(section_dict, required_opts)
    if not raise_error:
        return missing_opts
    elif len(missing_opts) > 0:
        msg = 'Section "{}" in the configuration is missing the following options: {}'.format(
            section, ', '.join(missing_opts)
        )
        raise ConfigurationError(msg)


def list_missing_subopts(required_subopts, config, section, option, raise_error=False):
    subopt_dict = config[section][option]
    missing_opts = _opts_missing_from_dict(subopt_dict, required_subopts)
    if not raise_error:
        return missing_opts
    elif len(missing_opts) > 0:
        msg = 'Option "{}" in section "{}" is missing the following sub-options: {}'.format(
            option, section, ', '.join(missing_opts)
        )
        raise ConfigurationError(msg)


def _opts_missing_from_dict(opt_dict, required_opts):
    missing_opts = []
    for opt in required_opts:
        if opt not in opt_dict.keys():
            missing_opts.append(opt)
    return missing_opts
