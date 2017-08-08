#!/usr/bin/env python3


#################
###### Need Implemented:
###### THERMAL_T2
###### unit tests
#################


"""
Generate mechanism solver file from a KPP-like mechanism file or one following PECANS style
"""

import argparse
from glob import glob
import os.path
import re

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

derivative_file = os.path.join(_mydir, 'chemderiv.pyx')
pyx_indent = '    '


############# SPECIFIC J Definitions ###############

MCM_J_Parameters = {
    #  in form L, M, N for use in MCM J equation :
    #  J = L * cosX ^ m * EXP(-N * secX), with X = solar zenith
    19 : [1.482e-06, 0.396, 0.298,1],
    14 : [2.879e-05, 1.067, 0.358,1],
    13 : [7.344e-06, 1.202, 0.417,1],
    20 : [7.600e-04, 0.396, 0.298,1],
    17 : [7.914e-05, 0.764, 0.364,1],
    7 : [2.644e-03, 0.261, 0.288,1],
    56 : [4.365e-05, 1.089, 0.323,1],
    54 : [4.095e-06, 1.111, 0.316,1],
    53 : [2.485e-06, 1.196, 0.328,1],
    55 : [1.135e-05, 0.974, 0.309]
    }

J_defins = {
    # R2smh Def : [TUV or MCM definition, which rate number, coefficient]
    'Pj_o33p': ['TUV_J', 3,1],
    'Pj_o31d': ['TUV_J',  2,1],
    'Pj_h2o2': ['TUV_J',  5,1],
    'Pj_no2': ['TUV_J',  6,1],
    'Pj_no3o2': ['TUV_J', 7,1],
    'Pj_no3o': ['TUV_J',  8,1],
    'Pj_hno3': ['TUV_J',  13,1],
    'Pj_hno4': ['TUV_J', 14,1],
    'Pj_ch2om': ['TUV_J',  19,1],
    'Pj_ch2or': ['TUV_J',  18,1],
    'Pj_glyc': ['TUV_J',  54,1],
    'Pj_bald': ['MCM_J', 19,1],
    'Pj_uald': ['MCM_J',  13,1],
    'Pj_ald': ['MCM_J',  14,1],
    'Pj_hpald': ['MCM_J', 20,1],
    'Pj_ibutald': ['MCM_J', 17,1],
    'Pj_hno2': ['MCM_J',  7,1],
    'Pj_noa': ['MCM_J', 56,1],
    'Pj_mpn': ['TUV_J', 27,1], # methyl peroxy nitrate
    'Pj_onit1': ['MCM_J', 53,1], # MCM n-propyl nitrate photolysis cross-section http://cohen.cchem.berkeley.edu/Portals/1/data/new_RONO2_mechanism.pdf
    'Pj_onitOH1': ['MCM_J', 53,(1/3)], # MCM n-propyl nitrate photolysis cross-section, divided by 3 due to the hydroxy group (Roberts and Fajer, 1989)
    'Pj_onitOH3': ['MCM_J', 55,(1/3)], # MCM tert-butyl nitrate photolysis cross-section, divided by 3 due to the hydroxy group (Roberts nd Fajer, 1989)
        # J = integral (quantum yeild * abs cross section * actinic flux dwavelength)
        # Jaltered = integral (quantum yeild/3 * abs cross section * actinic flux dwavelength)
        # which is 1/3 * integral (quantum yeild * abs cross section * actinic flux dwavelength)
        # and therefore Jaltered = J * 1/3. Nice and Simple!
    'Pj_iprno3': ['TUV_J', 33,1],
    'Pj_pana': ['TUV_J', 38,1],
    'Pj_panb': ['TUV_J', 39,1],
    'Pj_mvk': ['TUV_J', 44,1],
    'Pj_macr': ['TUV_J', 43,1],
    'Pj_ch3o2h': ['TUV_J', 24,1],
    'Pj_ch3cho': ['TUV_J', 20,1],
    'Pj_ch3cocho': ['TUV_J', 55,1],
    'Pj_ch3coo2h': ['TUV_J', 58,1],
    'Pj_ch3coch3': ['TUV_J', 48,1],
    'Pj_ch3coc2h5': ['TUV_J', 49,1],
    'Pj_hcochob': ['TUV_J', 54,1],
    'Pj_hcocho': ['TUV_J', 53,1],
    'Pj_hcochoc': ['TUV_J', 52,1]}

""" TUV photolysis Notes
   1 = O2 -> O + O
   2 = O3 -> O2 + O(1D)
   3 = O3 -> O2 + O(3P)
   4 = HO2 -> OH + O
   5 = H2O2 -> 2 OH
   6 = NO2 -> NO + O(3P)
   7 = NO3 -> NO + O2
   8 = NO3 -> NO2 + O(3P)
   9 = N2O -> N2 + O(1D)
  10 = N2O5 -> NO3 + NO + O(3P)
  11 = N2O5 -> NO3 + NO2
  12 = HNO2 -> OH + NO
  13 = HNO3 -> OH + NO2
  14 = HNO4 -> HO2 + NO2
  15 = NO3-(aq) -> NO2(aq) + O-
  16 = NO3-(aq) -> NO2-(aq) + O(3P)
  17 = NO3-(aq) with qy=1
  18 = CH2O -> H + HCO
  19 = CH2O -> H2 + CO
  20 = CH3CHO -> CH3 + HCO
  21 = CH3CHO -> CH4 + CO
  22 = CH3CHO -> CH3CO + H
  23 = C2H5CHO -> C2H5 + HCO
  24 = CH3OOH -> CH3O + OH
  25 = HOCH2OOH -> HOCH2O. + OH
  26 = CH3ONO2 -> CH3O + NO2
  27 = CH3(OONO2) -> CH3(OO) + NO2
  28 = CH3CH2ONO2 -> CH3CH2O + NO2
  29 = C2H5ONO2 -> C2H5O + NO2                           ethyl nitrate
  30 = n-C3H7ONO2 -> C3H7O + NO2                         n-prop nitrate
  31 = 1-C4H9ONO2 -> 1-C4H9O + NO2                       2-but nitrate
  32 = 2-C4H9ONO2 -> 2-C4H9O + NO2
  33 = CH3CHONO2CH3 -> CH3CHOCH3 + NO2
  34 = CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
  35 = CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
  36 = C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
  37 = C(CH3)3(ONO) -> C(CH3)3(O) + NO
  38 = CH3CO(OONO2) -> CH3CO(OO) + NO2
  39 = CH3CO(OONO2) -> CH3CO(O) + NO3
  40 = CH3CH2CO(OONO2) -> CH3CH2CO(OO) + NO2
  41 = CH3CH2CO(OONO2) -> CH3CH2CO(O) + NO3               ppn
  42 = CH2=CHCHO -> Products
  43 = CH2=C(CH3)CHO -> Products
  44 = CH3COCH=CH2 -> Products
  45 = HOCH2CHO -> CH2OH + HCO
  46 = HOCH2CHO -> CH3OH + CO
  47 = HOCH2CHO -> CH2CHO + OH
  48 = CH3COCH3 -> CH3CO + CH3
  49 = CH3COCH2CH3 -> CH3CO + CH2CH3
  50 = CH2(OH)COCH3 -> CH3CO + CH2(OH)
  51 = CH2(OH)COCH3 -> CH2(OH)CO + CH3
  52 = CHOCHO -> HCO + HCO
  53 = CHOCHO -> H2 + 2CO
  54 = CHOCHO -> CH2O + CO
  55 = CH3COCHO -> CH3CO + HCO
  56 = CH3COCOCH3 -> Products
  57 = CH3COOH -> CH3 + COOH
  58 = CH3CO(OOH) -> Products
  59 = CH3COCO(OH) -> Products
  60 = (CH3)2NNO -> Products  """

#####################################################

# Rates-related things
temperature_variable = 'TEMP'
ndens_air_variable = 'C_M'
time_variable = 'HOUR'
rate_expr_include_dir = os.path.join(_mydir, 'Rates')
c_math_fxns = ['exp', 'sqrt', 'log', 'log10', 'cos', 'abs']

inline_incl = False
inlines = ['cdef inline int int_min(int a, int b): return a if a <= b else b']

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

    def __init__(self, name, fixed=False):
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
        self.fixed = fixed

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

    @coefficient.setter
    def coefficient(self, value):
        self._coefficient = value

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
        if len(self.input_species) != 0:
            self._func_signature = 'd{}_dt(double {}, double {}, double {}, {})'.\
                format(self.changed_specie.name, time_variable, temperature_variable, ndens_air_variable,
                       ', '.join('double conc_'+s.name for s in self.input_species))
        else:
            self._func_signature = 'd{}_dt(double {}, double {})'.\
                format(self.changed_specie.name, temperature_variable, ndens_air_variable)

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

        if len(self.reactions) == 0:
            lines[-1] += '0'

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
            if not is_call and re.search('double\s+\w+', inpt) is None and re.search('int\s+\w+', inpt) is None:
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
            print(call_name)
            print(n_inputs)
            print(rate_ex.num_inputs)
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

def _parse_kpp_species(spec_file):
    """
    Parse a PECANS-style species file. The PECANS format requires that each specie is defined on
    its own line. Anything following a # is considered a comment and ignored
    :param spec_file: The path to the species file, as a str
    :return: nothing. All species are added to Species.instances.
    """
    with open(spec_file, 'r') as f:
        lineno = 0

        for line in f:
            lineno += 1

            if line.strip() == '#DEFVAR':
                fixed_species = False

            elif line.strip() == '#DEFFIX':
                fixed_species = True

            elif line[0] != '#':
                line = re.sub('\{(.*?)\}','', line) # Remove Comments
                if line.strip() == '': #Whole line was comment
                    pass
                else:
                    if len(re.findall('=', line)) < 1:
                        raise ReactionDefError('No reactant/product delimeter found on line {} of {}'.format(line, spec_file))
                    elif len(re.findall('=', line)) > 1:
                        raise ReactionDefError('Multiple reactant/product delimeters found on line {} of {}'.format(line, spec_file))

                    split_line = line.split('=')
                    spec_name = split_line[0]
                    spec_action = re.sub('\;', '', split_line[1])

                    if spec_action.strip() == 'IGNORE':
                        if len(spec_name) > 0:
                            Specie(spec_name.strip(), fixed_species) # don't need to return anything; automatically added to the "instances" list
                    else:
                           raise NotImplementedError('KPP action {} yet to be implmented'.format(spec_action.strip()))
            else:
                raise NotImplementedError('KPP Definition {} yet to be implmented'.format(line))

def _parse_kpp_reactions(rxn_file):
    """parse_pecan_reactions documentation:
    Parse a PECANS-style reaction file. The KPP format follows these rules:
    1) A '=' delineates the product and reactant sides of a reaction
    2) Within each side, individual species are separated by +
    3) Each species entry may include a coefficient before it; if omitted it is assumed to be one.
       The parser assumes that the species name begins at the first alphabetical character. It
       only allows numbers and decimal points (i.e. a period) in the coefficients
    4) The rate constant must follow a colon (:)
    5) Anything between {} is a comment and is ignored
        NOTE: // comments not implemented
    6) A semicolon designates the end of one line.
        NOTE: Mass Balance Checking is not implmented,
        only Chemical = IGNORE; definitions are currently accepted.
    7) Chemicals defined after DEFVAR are allowed to vary while those declared after DEFFIX
        are given a derivative of zero.
    :param rxn_file: the file to be parsed (as a string)
    :return: a list of reactions
    """

    reactions = []
    lineno = 0

    with open(rxn_file, 'r') as f:
        for line in f:
            lineno += 1
            line = re.sub('\{(.*?)\}','', line) # remove comments contained within {}

            if line.strip() == '#EQUATIONS':
                pass

            elif line[0] != '#':

                m = re.findall('=', line)

                if len(m) < 1:
                    raise ReactionDefError('No reactant/product delimeter found on line {} of {}'.format(lineno, rxn_file))

                elif len(m) > 1:
                    raise ReactionDefError('Multiple reactant/product delimeters found on line {} of {}'.format(lineno, rxn_file))

                i_rp = line.index(m[0])

                try:
                    i_rate = line.index(':')

                except ValueError:
                    raise ReactionDefError('No rate constant separator, :, found on line {} of {}'.format(lineno, rxn_file))

                react_str = re.sub('\+hv','',line[:i_rp].strip())
                prod_str = line[i_rp+1:i_rate].strip()
                rate_str = re.sub('\;','',line[i_rate+1:].strip())

                try:
                    reactants = [spec for spec in _iter_react_prod(react_str)]
                    products = [spec for spec in _iter_react_prod(prod_str)]
                except ReactionDefError as err:
                    raise ReactionDefError('Problem parsing reactants/products of line {} of {}: {}'.format(lineno, rxn_file, err.args[0])) from None


                rate_str = re.sub('_dp','',rate_str)
                rate_str = re.sub('D','e',rate_str)
                rate_str = re.sub('EXP','exp',rate_str)

                def self_evaluating(expr):
                    C_M, TEMP = 2.5e19, 300
                    try:
                        float(eval(expr))
                        return True
                    except NameError:
                        return False

                def nested_rates(rate_str):
                    # If it's a coefficient it can bet left alone
                    # if it is self evaluating
                    if not self_evaluating:
                        rate_str = rate_str[1:-1]
                        for rate_fraction in rate_str.split('*'):
                            if rate_str.startswith('('):
                                nested_rates(rate_fraction)
                            elif not self_evaluating(rate_fraction) and not rate_str.startswith('j'):
                                try:
                                    RateExpression.mark_rate_as_needed(rate_fraction)
                                except RateDefError as err:
                                    raise ReactionDefError('Problem parsing reactants/products of line {} of {}: {}'
                                                           .format(lineno, rxn_file, err.args[0])) from None

                if rate_str.startswith('j') or ('*' in rate_str and (rate_str.split('*')[1]).startswith('j')):
                    rate_str = define_photolysis(rate_str)

                if rate_str.startswith('('):
                    nested_rates(rate_str)

                else:
                    # If it's just a number that can be recognized as a float, leave it how it is
                        if not self_evaluating(rate_str):
                            try:
                                RateExpression.mark_rate_as_needed(rate_str)
                            except RateDefError as err:
                                raise ReactionDefError('Problem parsing reactants/products of line {} of {}: {}'
                                                       .format(lineno, rxn_file, err.args[0])) from None

                reactions.append(Reaction(reactants, products, rate_str))

            else:
                raise NotImplementedError('Unexpected Argument'.format(line))

    return reactions

def define_photolysis(rate_str):

    if '*' in rate_str:
        j_rate = rate_str.split('*')[1]
        coef = rate_str.split('*')[0]
    else:
        j_rate = rate_str
        coef = None

    j_arg = str(re.sub('\Aj','',j_rate))
    j_arg = str(re.sub('[\(\)]','',j_arg)).strip()

    j_func = J_defins[j_arg][0]
    n_rxn = J_defins[j_arg][1]
    coef2 = J_defins[j_arg][2]

    if not any([j_arg == key for key in list(J_defins.keys())]):
        raise RateDefError('Unrecognized Photolysis expression arguments {}'.format(j_arg))

    if j_func == 'MCM_J':
        args = ','.join([str(element) for element in MCM_J_Parameters[n_rxn]])
        func = 'MCM_J'
        if coef2 == 1:
            j_expression = '{}({},{})'.format(j_func, args, time_variable)
        else:
            j_expression = '{} * {}({},{})'.format(coef2, j_func, args, time_variable)
    elif j_func == 'TUV_J':
        if coef2 == 1:
            j_expression = '{}({},{})'.format(j_func, n_rxn, time_variable)
        else:
            j_expression = '{} * ({},{})'.format(coef2, j_func, n_rxn, time_variable)


    if coef:
        return '{}*{}'.format(coef, j_expression)
    else:
        return j_expression

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

def generate_kpp_mechanism(species_file, reactions_file, extra_rate_def_file=None):
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
    print('Building p')

    _read_rate_def_files(extra_rate_def_file)
    _parse_kpp_species(species_file)
    return _parse_kpp_reactions(reactions_file)

def generate_chemderiv_file(reactions):
    """
    Generates the chemical derivative .pyx file that contains all the derivative functions
     and rate constant functions plus the interface function that should be called from
     another Python program.
    :param reactions: a list of instances of Reactions
    :return: nothing
    """
    derivs = []
    for spec in Specie.instances:
        dC = Derivative(spec)
        if spec.fixed == False:
            for rxn in reactions:
                dC.add_rxn_if_relevant(rxn)
        dC.finalize()
        derivs.append(dC)

    with open(derivative_file, 'w') as f:
        f.write('from libc.math cimport {}\n\n'.format(', '.join(c_math_fxns)))
        if inline_incl:
            f.write('\n'.join([inline for inline in inlines]))
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

    lines.append('def chem_solver(double dt, double {}, double {}, double {}, {}):'.
        format(time_variable, temperature_variable, ndens_air_variable,
        ', '.join(['double conc_{}'.format(s.name) for s in Specie.instances]
                  )))

    for d in derivatives:
        lines.append(pyx_indent + 'cdef double dconc_{} = '.format(d.changed_specie.name) +
                     d.func_signature.replace('double', '') + ' * dt')
        dict_str.append('"conc_{0}": conc_{0} + dconc_{0}'.format(d.changed_specie.name))

    lines.append(pyx_indent + 'return {{ {0} }}'.format(', '.join(dict_str)))
    return lines

# Add styles and their respective parsing functions to this dictionary so that there
# is a single object describing which parsing function to call for which style

style_parse_fxn_dict = {'pecans':generate_pecans_mechanism, 'kpp':generate_kpp_mechanism}

def generate_mechanism(mechanism_style, species_file, reactions_file, additional_rates_file=None):
    try:
        build_fxn = style_parse_fxn_dict[mechanism_style]
    except KeyError:
        raise NotImplementedError('The mechgen.style_parse_fxn_dict does not have a parsing function listed for style '
                                  '{}'.format(mechanism_style))

    reactions = build_fxn(species_file, reactions_file, additional_rates_file)
    generate_chemderiv_file(reactions=reactions)

def _main():
    """
    Main driver function for this program
    :return: nothing
    """
    args = _get_args()
    generate_mechanism(mechanism_style=args.style.lower(), species_file=args.species_file,
                       reactions_file=args.reactions_file)

if __name__ == '__main__':
    _main()
