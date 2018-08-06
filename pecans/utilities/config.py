from ast import literal_eval
import configparser
import copy
import os
import re

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


class BetterConfig(configparser.RawConfigParser):
    """
    Extension of :class:`configparser.RawConfigParser` customized for use with PECANS

    This class adds a :func:`~pecans.utilities.config.BetterConfig.section_as_dict()` method that allows you to quickly
    return a section of the parser as a dictionary. It also overrides the :func:`configparser.RawConfigParser.read`
    method to parse input values from strings into Python types using the function set as the attribute `valuexform`.
    """
    # Based off of RawConfigParser to allow storing values not as strings
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.valuexform = str

    def section_as_dict(self, section):
        """
        Return the requested section's options as a dictionary
        :param section: the section name
        :type section: str

        :return: the dictionary representation of the section
        :rtype: dict
        """
        return {k: self[section][k] for k in self[section].keys()}

    def read(self, filenames, encoding=None):
        """
        Read one or more configuration files and stores their values internally
        :param filenames: the file or files to read
        :type filenames: str or list
        :param encoding: any encoding recognized by :func:`open`
        :return: the list of files successfully read
        """
        rvalue = super().read(filenames, encoding=encoding)
        for section in self.keys():
            for k in self[section].keys():
                self.set(section, k, self.valuexform(self.get(section, k)))

        return rvalue

    def __copy__(self):
        new_copy = self.__class__()
        for section in self.sections():
            new_copy.add_section(section)
            for opt in self.options(section):
                new_copy.set(section, opt, self.get(section, opt))
        return new_copy

    def __deepcopy__(self, memodict={}):
        new_copy = self.__class__()
        for section in self.sections():
            new_copy.add_section(section)
            for opt in self.options(section):
                new_copy.set(section, opt, copy.deepcopy(self.get(section, opt), memodict))
        return new_copy

    def as_string(self):
        s = ''
        for section in self.sections():
            s += '[{}]\n'.format(section)
            for opt in self.options(section):
                s += '{} = {}\n'.format(opt, self.get(section, opt))
            s += '\n'

        return s


def load_config_file(config_file):
    """
    Reads a configuration file and returns a BetterConfig instance representing it

    This function automatically instantiates a :class:`~pecans.utilities.config.BetterConfig` instance, sets it to parse
    values into Python literals (including tuples or dicts)
    :param config_file: the file to read.
    :type config_file: str or list
    :return: the :class:`~pecans.utilities.config.BetterConfig` object
    """
    if not os.path.isfile(config_file):
        raise IOError('Specified configuration file "{}" does not exist'.format(config_file))
    config = BetterConfig()
    config.valuexform = _parse_line_value
    config.read(config_file)
    return config


def _parse_as_dict(input_string):
    """
    Parse a value string from a config file into a dictionary

    This function takes a string of the form "key1: value1, key2: value2" etc. and parses it into a dictionary {key1:
    value1, key2: value2}. The string *must* have the key and value separated by a colon and key-value pairs separated
    by commas. If colons or commas are used in any other context, it will break. Individual values are parsed using
    :func:'~pecans.utilities.config._parse_as_literal`, keys are retained as strings.

    :param input_string: the string to parse
    :type input_string: str

    :return: the resulting dictionary
    :rtype: dict
    """
    # Parses an entry formatted as key : value, key : value, ...
    result = dict()
    entries = split_on_comma_re.split(input_string)
    for e in entries:
        try:
            key, val = e.split(':')
        except ValueError:
            raise ConfigurationError('Error parsing configuration file:\n'
                                     '  Attempting to parse it as a dictionary because of the presence of a colon,\n'
                                     '  but I found a piece that does not have exactly one colon:\n'
                                     '      {piece}'.format(piece=e))

        key = key.strip()
        val = _parse_as_literal(val.strip())

        result[key] = val

    return result


def _parse_as_tuple(input_string):
    """
    Parse a value string from a config file into a tuple

    This function takes a string of the form value1, value2, value3 and transforms it into a tuple: (value1, value2,
    value3). The values must be separated by commas; if commas are used in any other context, it cannot distinguish
    that. Individual values are parsed using :func:'~pecans.utilities.config._parse_as_literal`.

    :param input_string: the string to parse
    :type input_string: str

    :return: the resulting tuple
    :rtype: tuple
    """
    entries = split_on_comma_re.split(input_string)
    return tuple([_parse_as_literal(val.strip()) for val in entries])


def _parse_as_literal(val):
    """
    Try to parse the given string into a Python literal. The string must be exactly as it would be entered in the Python
    interpreter to be parsed. If it can't do so, it keeps it as a string.

    :param val: the string to parse
    :type val: str

    :return: the Python literal.
    """
    try:
        # literal_eval is a safer way of converting strings into Python objects, it will not evaluate arbitrary
        # expressions, only Python literal values. We can't use this to parse the whole entry because it won't be
        # formatted exactly as a Python dict.
        val = literal_eval(val.strip())
    except (ValueError, SyntaxError):
        # If it fails, we just keep it as a string
        pass

    return val


def _parse_line_value(input_string):
    """
    Overall parsing function that, given a value string from a config file, decides how to parse it.

    This function will hand off the input string to :func:'~pecans.utilities.config._parse_as_literal`,
    :func:'~pecans.utilities.config._parse_as_tuple`, or :func:'~pecans.utilities.config._parse_as_dict` depending on
    the format of the string (if it contains colons, it gets parsed as a dict; if commas but no colons, then as a tuple;
    otherwise just a regular Python literal value).

    :param input_string: the string to parse
    :type input_string: str

    :return: the parsed value
    """
    if ':' in input_string:
        return _parse_as_dict(input_string)
    elif ',' in input_string:
        return _parse_as_tuple(input_string)
    else:
        return _parse_as_literal(input_string)


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
    nx = config.get('DOMAIN', 'nx')
    ny = config.get('DOMAIN', 'ny')
    nz = config.get('DOMAIN', 'nz')

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
    section_dict = config.section_as_dict(section)
    missing_opts = _opts_missing_from_dict(section_dict, required_opts)
    if not raise_error:
        return missing_opts
    elif len(missing_opts) > 0:
        msg = 'Section "{}" in the configuration is missing the following options: {}'.format(
            section, ', '.join(missing_opts)
        )
        raise ConfigurationError(msg)


def list_missing_subopts(required_subopts, config, section, option, raise_error=False):
    subopt_dict = config.get(section, option)
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
