"""
Generate mechanism solver file from a KPP-like mechanism file or one following PECANS style
"""

import argparse
import re

derivative_file = 'chemderiv.pyx'
pyx_indent = '    '

class ChemError(Exception):
    pass

class SpeciesDefError(ChemError):
    pass

class ReactionDefError(ChemError):
    pass

class Specie:
    instances = []
    next_id = 0

    @property
    def name(self):
        return self._name

    @property
    def spec_id(self):
        return self._spec_id

    def __init__(self, name):
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
        cls.instances = []
        cls.next_id = 0

    @classmethod
    def find_by_name(cls, name, case_sensitive=False):
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
    next_id = 0

    @property
    def id(self):
        return self._id

    @property
    def reactants(self):
        return self._reactants

    @property
    def reactant_species(self):
        return [s.specie for s in self.reactants]

    @property
    def product_species(self):
        return [s.specie for s in self.products]

    @property
    def products(self):
        return self._products

    @property
    def rate_str(self):
        if isinstance(self._rate_fxn, float):
            return self._rate_fxn.__str__()
        else:
            return self._rate_fxn.__name__

    def __init__(self, reactants, products, rate_fxn):
        if not isinstance(reactants, (list, tuple, set)):
            raise TypeError('reactants must be an instance of list, tuple, or set')
        if not isinstance(products, (list, tuple, set)):
            raise TypeError('products must be an instance of list, tuple, or set')

        self._reactants = self.make_specie_tuple(reactants)
        self._products = self.make_specie_tuple(products)
        self._rate_fxn = rate_fxn
        self._id = self.__class__.next_id
        self.__class__.next_id += 1

    def get_reactant_specie(self, specie):
        for s in self.reactants:
            if s.specie == specie:
                return s
        return None

    def get_product_specie(self, specie):
        for s in self.products:
            if s.specie == specie:
                return s
        return None

    @classmethod
    def reset(cls):
        cls.instances = []
        cls.next_id = 0

    @staticmethod
    def make_specie_tuple(species):
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
    @property
    def changed_specie(self):
        return self._changed_specie

    @property
    def input_species(self):
        return self._input_species

    @property
    def reactions(self):
        return self._reactions

    @property
    def coefficients(self):
        return self._coefficients

    @property
    def func_signature(self):
        return self._func_signature

    def __init__(self, specie):
        if not isinstance(specie, Specie):
            raise TypeError('specie must be an instance of mechgen.Specie')

        self._changed_specie = specie
        self._input_species = set()
        self._reactions = []
        self._coefficients = []

    def add_rxn_if_relevant(self, rxn):
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
        # Transform the input species from a set to a tuple to ensure the ordering remains constant
        self._input_species = tuple(self._input_species)
        self._func_signature = 'd{}_dt({})'.format(self.changed_specie.name,
                                                   ', '.join('double conc_'+s.name for s in self.input_species))

    def func_def(self):
        if isinstance(self.input_species, set):
            raise RuntimeError('Derivative.func_def() called before Derivative.finalize()')

        lines = []
        lines.append('cdef double {}:'.format(self.func_signature))
        lines.append('    cdef double d{} = '.format(self.changed_specie.name))
        for i in range(len(self.reactions)):
            rxn = self.reactions[i]
            coeff = self.coefficients[i]
            rxn_str = '{}*{}*{}'.format(coeff, rxn.rate_str,
                                        '*'.join(['conc_'+s.name for s in rxn.reactant_species]))
            lines[-1] += rxn_str
            if i < len(self.reactions) - 1:
                lines[-1] += ' + '
        lines.append('    return d{}'.format(self.changed_specie.name))

        return lines


def get_args():
    parser = argparse.ArgumentParser(description='generates the needed .pyx files from a mechanism file')
    parser.add_argument('--solver', default='explicit', help='the numerical solver to use')
    parser.add_argument('--style', default='pecans', help='whether the input files are PECANS or KPP style')
    parser.add_argument('species_file', help='the file listing the species names')
    parser.add_argument('reactions_file', help='the file listing the chemical reactions and rate constants')

    args = parser.parse_args()
    return args


def parse_pecan_species(spec_file):
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


def parse_pecan_reactions(rxn_file):
    """parse_pecan_reactions documentation:
    parse_pecan_reactions(rxn_file)
        Parse a PECANS-style reaction file. The PECANS format follows these rules:
        1) Either -> or = delineates the product and reactant sides of a reaction
        2) Within each side, individual species are separated by +
        3) Each species entry may include a coefficient before it; if omitted it is assumed to be one.
           The parser assumes that the species name begins at the first alphabetical character. It
           only allows numbers and decimal points (i.e. a period) in the coefficients
        4) The rate constant must follow a colon (:)
        5) Anything after # is a comment and is ignored
        :param rxn_file: the file to be parsed (as a string)
        :return: nothing
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

            reactions.append(Reaction(reactants, products, 1e-15))

    return reactions


def _iter_react_prod(line):
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


def generate_chemderiv_file(reactions):
    derivs = []
    for spec in Specie.instances:
        dC = Derivative(spec)
        for rxn in reactions:
            dC.add_rxn_if_relevant(rxn)
        dC.finalize()
        derivs.append(dC)

    with open(derivative_file, 'w') as f:
        f.write('\n'.join(generate_interface_function(derivs)))
        f.write('\n\n')
        for d in derivs:
            f.write('\n'.join(d.func_def()))
            f.write('\n\n')


def generate_interface_function(derivatives):
    lines = []
    dict_str = []
    lines.append('def chem_solver(double dt, {}):'.format(', '.join(
        ['double conc_{}'.format(s.name) for s in Specie.instances]
    )))
    for d in derivatives:
        lines.append(pyx_indent + 'cdef double dconc_{} = '.format(d.changed_specie.name) +
                     d.func_signature.replace('double', '') + ' * dt')
        dict_str.append('"conc_{0}": conc_{0} + dconc_{0}'.format(d.changed_specie.name))

    lines.append(pyx_indent + 'return {{ {0} }}'.format(', '.join(dict_str)))
    return lines


def main():
    args = get_args()


if __name__ == '__main__':
    main()
