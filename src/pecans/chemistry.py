from collections.abc import Callable
from typing import TypeAlias
import numpy as np

from .ambient import AmbientConditions

RateFxn: TypeAlias = Callable[[AmbientConditions], float]

class Mechanism:
    """A representation of a chemical mechanism in memory

    An instance of this class is the end result of parsing mechanism files
    (or constructing a mechanism manually). Any mechanism file format will
    be represented by this class.

    Parameters
    ----------

    rate_functions
        A list of functions that each take a single argument that is a :class:`AmbientConditions` dataclass
        and return the rate or photolysis constant for the corresponding reaction given those ambient conditions.
        This list _must_ have a number of elements equal to the number of rows of ``reactant_coefficients``
        and be in the same order.
        
    reaction_coefficients
        An array that will have a number of rows equal to the number of reactions and a number of columns
        equal to the number of species in the mechanism. The values must indicate how many molecules of each
        species are consumed or produced by each reaction (i.e., their stoichiometric coefficients). Reactant
        consumed must have values < 0, products generated must have coefficients > 0, and species that do not
        participate in that reaction must have a value of 0.
    """
    def __init__(self, rate_functions: list[RateFxn], reaction_coefficients: np.ndarray) -> None:
        nrxn = len(rate_functions)
        if nrxn != reaction_coefficients.shape[0]:
            raise ValueError('There must be one rate function for each row of the reactant coefficient array')
        
        self.rate_functions = rate_functions
        self.ln_rate_constants = np.zeros(nrxn)
        self.reaction_coefficient_array = reaction_coefficients

        # Only the reactants are involved in the rate calculation. To generate the correct
        # matrix, we negate the values (so that reactants now have positive coefficients) and
        # then set any negative values to 0 (so that products do not contribute to the calculation).
        self.rate_coefficient_array = reaction_coefficients * -1
        self.rate_coefficient_array[self.rate_coefficient_array < 0] = 0.0

        # To avoid requiring ambient conditions during initialization, we tell the class that
        # it must calculate the rate constants at least the first time ``compute_rates`` is called.
        self.update_rate_constants = True

    def compute_rates(self, concentrations: np.ndarray, ambient_conditions: AmbientConditions, ambient_updated: bool) -> np.ndarray:
        # If the ambient conditions have changed or something else in the model invalidated the
        # rate constants, we must recalculated them. Usually self.update_rate_constants is only
        # true when the mechanism was just initialized and has not had this method called yet.
        if self.update_rate_constants or ambient_updated:
            for i, fxn in enumerate(self.rate_functions):
                k = fxn(ambient_conditions)
                self.ln_rate_constants[i] = np.log(k)
            self.update_rate_constants = False

        # The concentrations must be a column vector so that the dot product with the
        # rate coefficient array returns a column vector with n_reaction elements
        ln_concentration = np.log(concentrations.reshape(-1, 1))
        ln_rates = self.ln_rate_constants + self.rate_coefficient_array @ ln_concentration
        return np.exp(ln_rates)
    
    def compute_concentration_change(self, concentrations: np.ndarray, ambient_conditions: AmbientConditions, ambient_updated: bool) -> np.ndarray:
        rates = self.compute_rates(concentrations=concentrations, ambient_conditions=ambient_conditions, ambient_updated=ambient_updated)
        # The rates must be a row vector so that the dot product with the n_reaction by n_specie coefficient array
        # produces a vector with n_specie elements, which we then make sure only has 1 dimension to simplify the return value.
        rates = rates.reshape(1, -1)
        dc = rates @ self.reaction_coefficient_array
        return np.ravel(dc)
