from ast import literal_eval
import configparser

import pdb


class ConfigurationError(Exception):
    pass


class BetterConfig(configparser.RawConfigParser):
    # Based off of RawConfigParser to allow storing values not as strings
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.valuexform = str

    def section_as_dict(self, section):
        return {k: self[section][k] for k in self[section].keys()}

    def read(self, filenames, encoding=None):
        super().read(filenames, encoding=encoding)
        for section in self.keys():
            for k in self[section].keys():
                self.set(section, k, self.valuexform(self.get(section, k)))


def load_config_file(config_file):
    config = BetterConfig()
    config.valuexform = _parse_line_value
    config.read(config_file)
    return config


def _parse_as_dict(input_string):
    # Parses an entry formatted as key : value, key : value, ...
    result = dict()
    entries = input_string.split(',')
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
    entries = input_string.split(',')
    return tuple([_parse_as_literal(val.strip()) for val in entries])


def _parse_as_literal(val):
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
    if ':' in input_string:
        return _parse_as_dict(input_string)
    elif ',' in input_string:
        return _parse_as_tuple(input_string)
    else:
        return _parse_as_literal(input_string)


def get_domain_size_from_config(config):
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

    return shape