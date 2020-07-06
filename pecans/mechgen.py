#!/usr/bin/env python3

"""
Generate mechanism solver file from a KPP-like mechanism file or one following PECANS style
"""

#TODO: change so that the dict returned from chem solver is just the changes, not the new concentrations
#      (this should work better with op splitting)

import argparse
from glob import glob
import os.path
import re
from .utilities.config import load_config_file

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
config_file = os.path.join(os.path.dirname(__file__), '..', 'pecans_config.cfg')
config = load_config_file(config_file)
derivative_file = os.path.join(_mydir, 'chemderiv.pyx')
pyx_indent = '    '

# Rates-related things
temperature_variable = 'TEMP'
ndens_air_variable = 'CAIR'
rate_expr_include_dir = os.path.join(_mydir, 'Rates')
c_math_fxns = ['exp', 'sqrt', 'log', 'log10']


class ChemError(Exception):
    pass

class SpeciesDefError(ChemError):
    pass

class ReactionDefError(ChemError):
    pass

class RateDefError(ChemError):
    pass

class Specie:
    """
    Class which represents unique chemical species in the mechanism
    """
    instances = []
    next_id = 0

    @property
    def name(self):
        """
        The name of the species
        :return: the name as a string
        """
        return self._name

    @property
    def spec_id(self):
        """
        The numerical ID of the species
        :return: the ID as an integer
        """
        return self._spec_id

    def __init__(self, name):
        """
        Create a new, unique species. Each instance is automatically registered with the
        class variable "instances", and an error is thrown if the name matches one that
        already exists.
        :param name:
        :return: instance of Specie
        """
        if not isinstance(name, str):
            raise TypeError('name must be an instance of str')

        self._name = name
        self._spec_id = self.__class__.next_id
        other_names = [s.name.lower() for s in self.__class__.instances]
        if name.lower() in other_names:
            raise SpeciesDefError('name {} multiply defined'.format(name))

        self.__class__.instances.append(self)
        self.__class__.next_id += 1

    def __repr__(self):
        return '<{}.{} at {:#x}: {} (ID={})>'.format(__name__, self.__class__.__name__, id(self), self.name, self.spec_id)

    @classmethod
    def reset(cls):
        """
        Clears the instances list and resets the ID counter
        :return: nothing
        """
        cls.instances = []
        cls.next_id = 0

    @classmethod
    def find_by_name(cls, name, case_sensitive=False):
        """
        Finds the species instance by its name
        :param name: the name to search for
        :param case_sensitive: if true, matches case, if false (default), does not
        :return: an instance of Specie
        """
        if case_sensitive:
            for s in cls.instances:
                if s.name == name:
                    return s
        else:
            for s in cls.instances:
                if s.name.lower() == name.lower():
                    return s

        raise SpeciesDefError('Specie {} not defined'.format(name))


class ReactionSpecie:
    """
    Wrapper class around Specie that matches it with a coefficient
    """
    @property
    def specie(self):
        return self._specie

    @property
    def name(self):
        return self._specie.name

    @property
    def coefficient(self):
        return self._coefficient

    def __init__(self, specie, coefficient=1.0):
        """
        Create a new instance of ReactionSpecie that contains the given specie and
        its coefficient in the reaction.
        :param specie: An instance of Specie, representing the unique chemical specie in question
        :param coefficient: The coefficient (as a float). Defaults to 1.0 if not given
        :return: instance of ReactionSpecie
        """
        if not isinstance(specie, Specie):
            raise TypeError('specie must be of type mechgen.Specie')
        if not isinstance(coefficient, float):
            raise TypeError('coefficient must be of type float')

        self._specie = specie
        self._coefficient = coefficient

    def __str__(self):
        return '{} {}'.format(self.coefficient, self.name)

    def __repr__(self):
        return '<{}.{} at {:#x}: {}>'.format(__name__, self.__class__.__name__, id(self), self.__str__())


class Reaction:
    """
    Class representing chemical reactions in a mechanism
    """
    next_id = 0

    @property
    def id(self):
        """
        The unique numerical ID of the reaction
        :return: the ID as an integer
        """
        return self._id

    @property
    def reactants(self):
        """
        Returns the reactants of the reaction
        :return: a tuple of instances of ReactionSpecie
        """
        return self._reactants

    @property
    def reactant_species(self):
        """
        Returns the unique species present as reactants
        :return: a tuple of instances of Specie
        """
        return tuple([s.specie for s in self.reactants])

    @property
    def products(self):
        """
        Returns the products of the reaction
        :return: a tuple of instances of ReactionSpecie
        """
        return self._products

    @property
    def product_species(self):
        """
        Returns the unique species present as products
        :return: a tuple of instances of Specie
        """
        return tuple([s.specie for s in self.products])

    @property
    def rate_str(self):
        """
        Returns a representation of the rate constant as a string
        NOTE: behavior may change
        :return: a string, a number as a string or the name of the function given
        """
        if isinstance(self._rate_fxn, float):
            return self._rate_fxn.__str__()
        else:
            return self._rate_fxn

    def __init__(self, reactants, products, rate_fxn):
        """
        Create an instance of Reaction that contains the reactants, products, and rate constant
        :param reactants: a list, tuple, or set of instances of ReactionSpecie representing the reactants
        :param products: a list, tuple, or set of instances of ReactionSpecie representing the products
        :param rate_fxn: a representation of the rate constant, either as a function or value
        :return: instance of Reaction
        """
        if not isinstance(reactants, (list, tuple, set)):
            raise TypeError('reactants must be an instance of list, tuple, or set')
        elif any([not isinstance(reactant, ReactionSpecie) for reactant in reactants]):
            raise TypeError('All elements of reactants must be instances of ReactionSpecie')
        if not isinstance(products, (list, tuple, set)):
            raise TypeError('products must be an instance of list, tuple, or set')
        elif any([not isinstance(product, ReactionSpecie) for product in products]):
            raise TypeError('All elements of products must be instances of ReactionSpecie')

        self._reactants = self._make_specie_set(reactants)
        self._products = self._make_specie_set(products)
        self._rate_fxn = rate_fxn
        self._id = self.__class__.next_id
        self.__class__.next_id += 1

    def get_reactant_specie(self, specie):
        """
        Given an instance of Specie, returns the first corresponding reactant ReactionSpecie of that type of Specie
        :param specie: an instance of Specie representing the chemical specie in question
        :return: an instance of ReactionSpecie
        """
        for s in self.reactants:
            if s.specie == specie:
                return s
        return None

    def get_product_specie(self, specie):
        """
        Given an instance of Specie, returns the first corresponding product ReactionSpecie of that type of Specie
        :param specie: an instance of Specie representing the chemical specie in question
        :return: an instance of ReactionSpecie
        """
        for s in self.products:
            if s.specie == specie:
                return s
        return None

    @classmethod
    def reset(cls):
        """
        Restores the class variables to their initial state
        :return: nothing
        """
        cls.next_id = 0

    @staticmethod
    def _make_specie_set(species):
        """
        (internal) convert a list, tuple, or set of instances of ReactionSpecie to a set with no repeated specie types
        If any species are repeated multiple times in the input, their coefficients are added together
        to give a single representation in the output. If given a set, does nothing.
        :param species: a list, tuple, or set of instances of ReactionSpecie
        :return: a set of instances ReactionSpecie
        """
        if isinstance(species, set):
            return species

        spec_dict = {}
        for s in species:
            if s.name not in spec_dict.keys():
                spec_dict[s.name] = s
            else:
                spec_dict[s.name].coefficient += s.coefficient

        return set(spec_dict.values())

    def __str__(self):
        return '{} -> {} : {}'.format(
            ' + '.join(['{}*{}'.format(s.coefficient, s.name) for s in self.reactants]),
            ' + '.join(['{}*{}'.format(s.coefficient, s.name) for s in self.products]),
            self._rate_fxn
        )

    def __repr__(self):
        return '<{}.{} at {:#x}: {}>'.format(__name__, self.__class__.__name__, id(self), self.__str__())


class Derivative:
    """
    Represents the derivative needed to calculate the change in a chemical specie
    """
    max_summands_per_line = 20

    @property
    def changed_specie(self):
        """
        The specie that the derivative is for
        :return: an instance of Specie
        """
        return self._changed_specie

    @property
    def input_species(self):
        """
        The species whose concentrations are needed as inputs to calculate the derivative
        :return: a set of instances of Specie
        """
        return self._input_species

    @property
    def reactions(self):
        """
        The reactions that involve the changed_specie
        :return: a list of instances of Reaction
        """
        return self._reactions

    @property
    def coefficients(self):
        """
        The coefficients for each reaction, tracking both how many molecules of the changed_specie
        are produced or consumed and the sign (positive for produced, negative for consumed). This
        is an overall coefficient, so if a specie appears on both sides, there will still be only
        one coefficient for that reaction.
        :return: a list of floats
        """
        return self._coefficients

    @property
    def func_signature(self):
        """
        The signature of the function that will be called to calculate the derivative. Includes the
        function name, argument types, and argument names, i.e. dNO_dt(double conc_NO, double conc_O3)
        Only defined after calling finalize() on this instance.
        :return: a string that represents the function signature
        """
        return self._func_signature

    def __init__(self, specie):
        """
        Initialize an instance of Derivative that can be used to represent the change of a chemical
        specie in a single timestep.
        :param specie: The specie that the derivative is for; an instance of Specie
        :return: an instance of Derivative
        """
        if not isinstance(specie, Specie):
            raise TypeError('specie must be an instance of mechgen.Specie')

        self._changed_specie = specie
        self._input_species = set()
        self._reactions = []
        self._coefficients = []

    def add_rxn_if_relevant(self, rxn):
        """
        Given a reaction, this will determine if the changed_specie appears on either side of the reaction
        and if so, adds the reaction to the list of information that the derivative will need
        :param rxn: an instance of Reaction
        :return: a boolean, True if the reaction was relevant, False otherwise
        """
        did_add = False
        if isinstance(self.input_species, tuple):
            raise RuntimeError('This instance of Derivative has already been finalized and cannot be modified further')

        as_prod = rxn.get_product_specie(self.changed_specie)
        as_react = rxn.get_reactant_specie(self.changed_specie)
        if as_prod is None:
            prod_coeff = 0
        else:
            prod_coeff = as_prod.coefficient

        if as_react is None:
            react_coeff = 0
        else:
            react_coeff = as_react.coefficient

        if as_react is not None or as_prod is not None:
            did_add = True
            self._input_species |= set(rxn.reactant_species)
            self._reactions.append(rxn)
            self._coefficients.append(prod_coeff - react_coeff)

        return did_add

    def finalize(self):
        """
        To be called once all reactions are added. This, internally, converts the collection of input species
        from a set to a tuple, to ensure that the order is consistent at all times.
        :return: nothing
        """
        # Transform the input species from a set to a tuple to ensure the ordering remains constant
        self._input_species = tuple(self._input_species)
        self._func_signature = 'd{}_dt(double {}, double {}, {})'.\
            format(self.changed_specie.name, temperature_variable, ndens_air_variable,
                   ', '.join('double conc_'+s.name for s in self.input_species))

    def func_def(self):
        """
        Creates the Cython function definition for the function that will calculate this derivative.
        :return: A list of strings, each string is a line in the function definition
        """
        if isinstance(self.input_species, set):
            raise RuntimeError('Derivative.func_def() called before Derivative.finalize()')

        lines = []
        # Create the function signature
        lines.append('cdef double {}:'.format(self.func_signature))
        # The next line defines the return variable
        return_var = 'd{}'.format(self.changed_specie.name)
        lines.append('{}cdef double {} = '.format(pyx_indent, return_var))
        for i in range(len(self.reactions)):
            # We then loop through the reactions that change this specie and sum each one up
            # A quirk of how Cython seems to handle addition is that it breaks a line down into
            # operations and recursively assigns summands as objects, so a long equation can have
            # summand objects that are too deeply nested. We get around this by separating long
            # expressions over multiple lines
            if i % Derivative.max_summands_per_line == 0 and i > 0:
                lines.append('{}{} += '.format(pyx_indent, return_var))

            rxn = self.reactions[i]
            coeff = self.coefficients[i]
            rxn_str = '{}*{}*{}'.format(
                coeff, rxn.rate_str, '*'.join(
                    ['conc_'+s.name if s.coefficient == 1.0 else 'conc_{}**{}'.format(
                        s.name, s.coefficient
                    ) for s in rxn.reactants]
                )
            )
            lines[-1] += rxn_str
            if i < len(self.reactions) - 1 and (i+1) % Derivative.max_summands_per_line:
                # There are only two cases where we do not want to add a "+" sign after each reaction's
                # contribution: if it's the last reaction, or if we are going to split the summation onto
                # the next line
                lines[-1] += ' + '
        lines.append('    return d{}'.format(self.changed_specie.name))

        return lines


class RateExpression:
    """
    Class that holds all the defined rate constant expression functions
    """
    instances = []

    def __init__(self, rate_string, rate_file, file_line_number):
        """
        Create a new instance of RateExpression that holds information about a rate expression function.
        That instance is stored in the class variable instances (RateExpression.instances), so there is
        no need to store the returned instance.
        :param rate_string: the string containing the entire function definition for the rate expression.
        :param rate_file: the file that the definition was read from. Used to print more helpful error messages.
        :param file_line_number: the line number in the file where the definition began. Again, for helpful error
        messages.
        :return: instance of RateExpression, which is also stored in RateExpression.instances.
        """
        self.used_in_mech = False
        self.declaring_file = os.path.abspath(rate_file)
        self.file_line = file_line_number
        self.rate_name, self.num_inputs, self.temperature_ind, self.ndens_air_ind = RateExpression._check_rate_string(rate_string)
        self.body = rate_string
        for rate_ex in RateExpression.instances:
            if self.rate_name == rate_ex.rate_name:
                raise RateDefError('Rate expression named "{}" already exists at line {} of file {}'.
                                   format(self.rate_name, rate_ex.file_line, rate_ex.declaring_file))
        RateExpression.instances.append(self)

    @classmethod
    def _check_rate_string(cls, rate_string, is_call=False):
        """
        Verify that the given rate string is a valid rate expression function.
        :param rate_string: the string containing at least the cdef <name>( args ): bit.
        :param is_call: optional boolean, default False. If False, this function assumes that the string is a definition
        of the function, and checks for things like the "cdef" keyword and "double" type for all inputs. If True, then
        it is assumed that the string is a call to the function rather than the definition, and things not relevant to a
        call are ignored. This allows one to use this function to parse either a definition or call and verify that the
        call has the same characteristics as the definition.
        :return: four values:
            function_name - the name of the function
            n_inputs - how many inputs to the function
            index_temperature - which of the inputs is the temperature value (0 based)
            index_air_number_density - which of the inputs is the number density of air (0 based)
        """
        function_def = rate_string.split(':')[0]

        m = re.match('[cp]*def', rate_string)
        if not is_call and (m is None or m.group() != 'cdef'):
            raise RateDefError('Rate expressions must be defined as cdef functions')

        m = re.search('\w+(?=\()', function_def)
        if m is None:
            raise RateDefError('Could not identify rate expression name for string:\n{}'.format(rate_string))
        function_name = m.group()

        if not is_call:
            # If this is a definition, not a call, then the inputs should be valid variable names, separated by
            # commas/spaces (meaning only letters, numbers, and underscores should be there, besides commas and spaces)
            m = re.search('(?<=\()[\w\s,]+(?=\))', function_def)
        else:
            # If it is a call in the equations file, then there could be any manner of symbols in there
            m = re.search('(?<=\().+(?=\))', function_def)
        if m is None:
            raise RateDefError('Could not identify rate expression inputs for string:\n{}'.format(rate_string))
        function_inputs = m.group().split(',')
        index_temperature = None
        index_air_number_density = None
        for inpt in function_inputs:
            if not is_call and re.search('double\s+\w+', inpt) is None:
                raise RateDefError('Input "{}" is not defined as type "double"'.format(inpt))
            elif temperature_variable in inpt:
                index_temperature = function_inputs.index(inpt)
            elif ndens_air_variable in inpt:
                index_air_number_density = function_inputs.index(inpt)

        return function_name, len(function_inputs), index_temperature, index_air_number_density

    @classmethod
    def find_rate_by_name(cls, rate_name):
        """
        Given a rate expression name, finds the instance that corresponds to it.
        :param rate_name: The name of the rate expression, as a string
        :return: the instance of RateExpression with that name. Raises a KeyError if one is not found.
        """
        for rate_ex in cls.instances:
            if rate_ex.rate_name == rate_name:
                return rate_ex
        raise KeyError('No rate named "{}" defined'.format(rate_name))

    @classmethod
    def mark_rate_as_needed(cls, rate_call):
        """
        Marks that the specified rate expression is used in the current mechanism, which causes it to be
        inlined in the chemderiv.pyx file.
        :param rate_call: the call to the rate expression, e.g. ARR2( 1.2e-12, 1310.0, TEMP )
        :return: none
        """
        rate_name = cls._check_rate_call(rate_call)
        cls.find_rate_by_name(rate_name).used_in_mech = True

    @classmethod
    def _check_rate_call(cls, call_string):
        """
        Checks the given rate expression call against the definition stored in RateExpression.instance
        :param call_string: the call to the rate expression, e.g. ARR2( 1.2e-12, 1310.0, TEMP )
        :return: the name of the call, in the example for call_string, "ARR2"
        """
        call_name, n_inputs, i_temp, i_cair = cls._check_rate_string(call_string, is_call=True)
        try:
            rate_ex = cls.find_rate_by_name(call_name)
        except KeyError as err:
            raise RateDefError(err.args[0]) from None

        if n_inputs != rate_ex.num_inputs:
            raise RateDefError('Rate call has a different number of arguments ({}) than its definition ({})'.
                               format(n_inputs, rate_ex.num_inputs))
        elif i_temp != rate_ex.temperature_ind:
            raise RateDefError('Rate call has temperature in a different place ({}) than its definition ({})'.
                               format(i_temp + 1, rate_ex.temperature_ind + 1))
        elif i_cair != rate_ex.ndens_air_ind:
            raise RateDefError('Rate call has the number density of air in a different place ({}) than its definition ({})'.
                               format(i_cair + 1, rate_ex.ndens_air_ind + 1))

        return call_name


    """
    @classmethod
    def _cythonize_rate_expression(cls, rate_string):
        function_def, function_body = rate_string.split(':')
        m = re.search('\w+(?=\()', function_def)
        if m is None:
            raise RateDefError('Could not identify rate expression name for string:\n{}'.format(rate_string))
        function_name = m.group()

        m = re.search('(?<=\().+(?=\))')
        if m is None:
            raise RateDefError('Could not identify rate expression inputs for string:\n{}'.format(rate_string))
        function_inputs = [s.strip() for s in m.group().split(',')]
        for rinput in function_inputs:
            if re.search('\W', rinput) is not None:
                raise RateDefError('Invalid character in rate expression input "{}"'.format(rinput))
    """


def _get_args():
    """
    Get the command line arguments
    :return: an instance of an argparse namespace
    """
    parser = argparse.ArgumentParser(description='generates the needed .pyx files from a mechanism file')
    parser.add_argument('--solver', default='explicit', help='the numerical solver to use')
    parser.add_argument('--style', default='pecans', help='whether the input files are PECANS or KPP style')
    parser.add_argument('species_file', help='the file listing the species names')
    parser.add_argument('reactions_file', help='the file listing the chemical reactions and rate constants')

    args = parser.parse_args()
    return args


def _parse_pecan_species(spec_file):
    """
    Parse a PECANS-style species file. The PECANS format requires that each specie is defined on
    its own line. Anything following a # is considered a comment and ignored
    :param spec_file: The path to the species file, as a str
    :return: nothing. All species are added to Species.instances.
    """
    with open(spec_file, 'r') as f:
        for line in f:
            try:
                i = line.index('#')
            except ValueError:
                line2 = line.strip()
            else:
                line2 = line[:i].strip()

            if len(line2) > 0:
                Specie(line2)  # don't need to return anything; automatically added to the "instances" list


def _parse_pecan_reactions(rxn_file):
    """parse_pecan_reactions documentation:
    Parse a PECANS-style reaction file. The PECANS format follows these rules:
    1) Either -> or = delineates the product and reactant sides of a reaction
    2) Within each side, individual species are separated by +
    3) Each species entry may include a coefficient before it; if omitted it is assumed to be one.
       The parser assumes that the species name begins at the first alphabetical character. It
       only allows numbers and decimal points (i.e. a period) in the coefficients
    4) The rate constant must follow a colon (:)
    5) Anything after # is a comment and is ignored
    :param rxn_file: the file to be parsed (as a string)
    :return: a list of reactions
    """
    lineno = 0
    reactions = []
    with open(rxn_file, 'r') as f:
        for line in f:
            lineno += 1
            try:
                i = line.index('#')
            except ValueError:
                line2 = line.strip()
            else:
                line2 = line[:i].strip()

            # Allow reactants and products to be separated by -> or =
            m = re.findall('->|=', line2)
            if len(m) < 1:
                raise ReactionDefError('No reactant/product delimeter found on line {} of {}'.format(lineno, rxn_file))
            elif len(m) > 1:
                raise ReactionDefError('Multiple reactant/product delimeters found on line {} of {}'.format(lineno, rxn_file))

            i_rp = line2.index(m[0])

            # Find the beginning of the rate coefficient
            try:
                i_rate = line2.index(':')
            except ValueError:
                raise ReactionDefError('No rate constant separator, :, found on line {} of {}'.format(lineno, rxn_file))

            react_str = line2[:i_rp].strip()
            prod_str = line2[i_rp+1:i_rate].strip()
            rate_str = line2[i_rate+1:].strip()

            try:
                reactants = [spec for spec in _iter_react_prod(react_str)]
                products = [spec for spec in _iter_react_prod(prod_str)]
            except ReactionDefError as err:
                # https://stackoverflow.com/questions/33809864/disable-exception-chaining-in-python-3
                raise ReactionDefError('Problem parsing reactants/products of line {} of {}: {}'
                                       .format(lineno, rxn_file, err.args[0])) from None

            if rate_str.startswith('j'):
                pass # deal with photolysis
            else:
                # If it's just a number that can be recognized as a float, leave it how it is
                try:
                    float(rate_str)
                except ValueError:
                    try:
                        RateExpression.mark_rate_as_needed(rate_str)
                    except RateDefError as err:
                        raise ReactionDefError('Problem parsing reactants/products of line {} of {}: {}'
                                               .format(lineno, rxn_file, err.args[0])) from None

            reactions.append(Reaction(reactants, products, rate_str))

    return reactions


def _iter_react_prod(line):
    """
    Iterates over reactants OR products, not both. Assumes that they are split by + signs
    :param line: The str of reactants or products read in from the reactions file
    :return: is an iterable
    """
    for spec in line.split('+'):
        m = re.search('[a-zA-Z]', spec)
        if m is None:
            raise ReactionDefError('Could not find a species name in substring "{}"'.format(spec))

        coeff_str = spec[:m.start()].strip()
        name = spec[m.start():].strip()

        if len(coeff_str) > 0:
            try:
                coeff = float(coeff_str)
            except ValueError:
                raise ReactionDefError('Could not interpret "{}" as a numeric value'.format(coeff))
        else:
            coeff = 1.0

        if coeff < 0.0:
            raise ReactionDefError('Negative coefficient detected')

        try:
            specie = Specie.find_by_name(name)
        except SpeciesDefError:
            raise ReactionDefError('Species name "{}" not defined'.format(name))

        yield ReactionSpecie(specie, coeff)


def _read_rate_def_files(additional_file):
    """
    Reads in all *.rate files in the Rates subdirectory, plus any files specified in the additional_file list
    :param additional_file: a list of additional files that have rate expression definitions
    :return: none. All rate expressions are added to RateExpression.instances.
    """
    if additional_file is not None:
        if not isinstance(additional_file, str):
            raise TypeError('additional_file must be a string')
        elif not os.path.isfile(additional_file):
            raise IOError('additiona_file ({}) does not exist'.format(additional_file))

    rate_files = glob(os.path.join(rate_expr_include_dir, '*.rate'))
    if additional_file is not None:
        rate_files.append(additional_file)

    for rfile in rate_files:
        if not os.path.isfile(rfile):
            raise IOError('Rates file {} does not exist'.format(rfile))

        reading_rate = False
        curr_rate_expr = []
        with open(rfile, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if re.search('def.+:', line) is not None:
                    reading_rate = True
                    if len(curr_rate_expr) > 0:
                        RateExpression(''.join(curr_rate_expr), rfile, def_line_num)
                    curr_rate_expr = [line]
                    def_line_num = line_num
                elif reading_rate and len(line.strip()) > 0:
                    curr_rate_expr.append(line)
            # Add the last rate expression in the file
            RateExpression(''.join(curr_rate_expr), rfile, def_line_num)


def generate_pecans_mechanism(species_file, reactions_file, extra_rate_def_file=None):
    """
    Main function that generates a mechanism file for the PECANS style inputs. Any other
    input file style must have an equivalent primary function that reads rate expression
    files, reads the species file, and reads the reactions file.
    :param species_file: the file defining all the species included in the mechanism
    :param reactions_file: the file defining the reactions that comprise the mechanism
    :param *extra_rate_def_file: any additional files beyond those in the "Rates" subdirectory
     that define rate constant expressions.
    :return:
    """
    _read_rate_def_files(extra_rate_def_file)
    _parse_pecan_species(species_file)
    return _parse_pecan_reactions(reactions_file)


def generate_chemderiv_file(reactions, additional_params):
    """
    Generates the chemical derivative .pyx file that contains all the derivative functions
     and rate constant functions plus the interface function that should be called from
     another Python program.
    :param reactions: a list of instances of Reactions
    :return: nothing
    :param ode_solver: the choice of chemistry kinetics solver
    :return: nothing
    :param additional_params: the additional parameters for chemderive, for instance, alpha value
    :return: nothing
    """
    derivs = []
    for spec in Specie.instances:
        dC = Derivative(spec)
        for rxn in reactions:
            dC.add_rxn_if_relevant(rxn)
        dC.finalize()
        derivs.append(dC)

    with open(derivative_file, 'w') as f:
        f.write('from libc.math cimport {}\n'.format(', '.join(c_math_fxns)))
        f.write('import numpy as np\n\n')
        f.write('{}'.format('\n'.join(additional_params)))
        f.write('\n\n')
        f.write('\n'.join(_generate_interface_ode_function(derivs)))
        f.write('\n\n')
        f.write('\n'.join(_generate_interface_function(derivs)))
        f.write('\n\n')
        for d in derivs:
            f.write('\n'.join(d.func_def()))
            f.write('\n\n\n')

        rate_ex_comment = '# RATE EXPRESSIONS #'
        f.write('{0}\n{1}\n{0}\n\n'.format('#'*len(rate_ex_comment), rate_ex_comment))
        for rate_ex in RateExpression.instances:
            if rate_ex.used_in_mech:
                f.write(rate_ex.body)
                f.write('\n\n')


def _generate_interface_function(derivatives):
    """
    Helper function that generates the interface function
    :param derivatives: a list of instances of Derivatives that represents all the derivatives that
    need to be calculated
    :return: a list of strings, each string is a line in the interface function
    """
    lines = []
    dict_str = []
    lines.append('def chem_solver(double dt, double {}, double {}, {}):'.
        format(temperature_variable, ndens_air_variable,
        ', '.join(['double conc_{}'.format(s.name) for s in Specie.instances]
                  )))
    for d in derivatives:
        lines.append(pyx_indent + 'cdef double dconc_{} = '.format(d.changed_specie.name) +
                     d.func_signature.replace('double', '') + ' * dt')
        dict_str.append('"{0}": conc_{0} + dconc_{0}'.format(d.changed_specie.name))

    lines.append(pyx_indent + 'return {{ {0} }}'.format(', '.join(dict_str)))
    return lines

def _generate_interface_ode_function(derivatives):
    """
    Helper function that generates the interface function
    :param derivatives: a list of instances of Derivatives that represents all the derivatives that
    need to be calculated
    :param ode_solver: the choice of ode solver defined in assimulo package
    :return: a list of strings, each string is a line in the interface function
    """
    if 'const_species' in config.section_as_dict('CHEMISTRY'):
        const_species = list(config.get('CHEMISTRY', 'const_species').keys())
    else:
        const_species = list()
    if 'forced_species' in config.section_as_dict('CHEMISTRY'):
        const_species.extend(list(config.get('CHEMISTRY', 'forced_species')))

    lines = []
    dict_str = []

    lines.append('def rhs(f, double t, param):')
    lines.append(pyx_indent + '[{}] = f'.format(','.join(['this_conc_{} '.format(s.name) for s
                                                          in Specie.instances if s.name not in const_species])))
    lines.append(pyx_indent + '[TEMP, CAIR, {}] = param'.format(','.join(['this_conc_{} '.format(specie)
                                                                           for specie in const_species])))
    for d in derivatives:
        if d.changed_specie.name not in const_species:
            lines.append(pyx_indent + 'cdef double dconc_{}_dt = '.format(d.changed_specie.name) +
                    d.func_signature.replace('double', '').replace('conc', 'this_conc'))
            dict_str.append('dconc_{0}_dt'.format(d.changed_specie.name))

    lines.append(pyx_indent + 'return np.array([{0}])'.format(', '.join(dict_str)))
    return lines

# Add styles and their respective parsing functions to this dictionary so that there
# is a single object describing which parsing function to call for which style

style_parse_fxn_dict = {'pecans':generate_pecans_mechanism}

def generate_mechanism(mechanism_style, species_file, reactions_file, additional_rates_file=None,
                       additional_params=None):
    try:
        build_fxn = style_parse_fxn_dict[mechanism_style]
    except KeyError:
        raise NotImplementedError('The mechgen.style_parse_fxn_dict does not have a parsing function listed for style '
                                  '{}'.format(mechanism_style))

    reactions = build_fxn(species_file, reactions_file, additional_rates_file)
    generate_chemderiv_file(reactions=reactions, additional_params=additional_params)

def _main():
    """
    Main driver function for this program
    :return: nothing
    """
    args = _get_args()
    generate_mechanism(mechanism_style=args.style.lower(), species_file=args.species_file,
                       reactions_file=args.reactions_file, ode_solver=args.ode_solver)


if __name__ == '__main__':
    _main()
