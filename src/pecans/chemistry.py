from collections.abc import Callable
from typing import TypeAlias
import numpy as np

from .ambient import AmbientConditions

RateFxn: TypeAlias = Callable[[AmbientConditions], float]

# TODO: to make these functions broadcast concentrations

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

    def compute_rates(self, concentrations: np.ndarray, ambient_conditions: AmbientConditions, ambient_updated: bool, squeeze: bool = True) -> np.ndarray:
        """Compute the rates of all reactions for all grid cells.

        The rates will give how many times that reaction occurs per unit time.

        Parameters
        ----------
        concentrations
            The array of concentrations for all species for all grid cells. It must have the different species
            as the final dimension. That is, for a 1D model with 10 grid cells and 4 species, this must be 10-by-4.
            For a 2D model with nx = 10 and ny = 5 (and 4 species again), this would be 10-by-5-by-4.

        ambient_conditions
            An instance of a class that provides the necessary ambient conditions (temperature, pressure, etc.) needed
            to calculate the rate constants. Note that standard kinetic rate constants typically need temperature and
            pressure, but photolysis rates may be parameterized a variety of ways. Thus, it is important that this
            type provide all of the conditions required by the rate functions passed in when the Mechanism instance
            was created.

        ambient_updated
            Set this to ``True`` if the ambient conditions have changed since the last time step, requiring the
            rate constants to be recalculated. If ``False``, the rate constants will not be recalculated. If your
            ambient conditions do change, the safest approach is to always pass ``True`` and accept the performance
            hit.

        squeeze
            If ``True`` (default), the returned array will have the same number of dimensions as ``concentrations``,
            but with the last dimension now equal in length to the number of reactions, rather than species. Setting
            this to ``False`` keeps an extra length-1 dimension on the end of the array that is a result of the
            way the matrix multiplication is handled internally. Most use outside of the PECANS internals should
            leave this as ``True``.

        Returns
        -------
        rates
            An array with the rates for all of the reactions in each grid cell. See ``squeeze`` for information about
            the shape.
        """
        # If the ambient conditions have changed or something else in the model invalidated the
        # rate constants, we must recalculated them. Usually self.update_rate_constants is only
        # true when the mechanism was just initialized and has not had this method called yet.
        if self.update_rate_constants or ambient_updated:
            for i, fxn in enumerate(self.rate_functions):
                k = fxn(ambient_conditions)
                self.ln_rate_constants[i] = np.log(k)
            self.update_rate_constants = False

        # The rules for broadcasting matrix multiplication is that the matrices to multiply
        # must be in the last two dimensions. So if concentrations was nx-by-nspecies, we
        # need to make it nx-by-nspecies-by-1 so that the multiplication broadcasts over the
        # nx grid cells.
        ln_concentration = np.log(concentrations[..., np.newaxis])
        ln_rates = self.ln_rate_constants + self.rate_coefficient_array @ ln_concentration

        # ln_rates will have an extra singleton dimension at the end
        if squeeze:
            return np.exp(ln_rates.squeeze())
        else:
            return np.exp(ln_rates)

    def compute_concentration_change(self, concentrations: np.ndarray, ambient_conditions: AmbientConditions, ambient_updated: bool) -> np.ndarray:
        """Compute the change in concentrations of each specie in each grid cell.
        
        Parameters
        ----------
        concentrations
            The array of concentrations for all species for all grid cells. It must have the different species
            as the final dimension. That is, for a 1D model with 10 grid cells and 4 species, this must be 10-by-4.
            For a 2D model with nx = 10 and ny = 5 (and 4 species again), this would be 10-by-5-by-4.

        ambient_conditions
            An instance of a class that provides the necessary ambient conditions (temperature, pressure, etc.) needed
            to calculate the rate constants. Note that standard kinetic rate constants typically need temperature and
            pressure, but photolysis rates may be parameterized a variety of ways. Thus, it is important that this
            type provide all of the conditions required by the rate functions passed in when the Mechanism instance
            was created.

        ambient_updated
            Set this to ``True`` if the ambient conditions have changed since the last time step, requiring the
            rate constants to be recalculated. If ``False``, the rate constants will not be recalculated. If your
            ambient conditions do change, the safest approach is to always pass ``True`` and accept the performance
            hit.

        Returns
        -------
        delta_conc
            The changes in concentrations, the same shape as ``concentrations``.
        """
        rates = self.compute_rates(concentrations=concentrations, ambient_conditions=ambient_conditions, ambient_updated=ambient_updated, squeeze=False)

        # If concentrations was nx-by-nspecies, then rates will now be nx-by-nreactions-by-1. We need to transpose the last two dimensions so that
        # we properly match the reactions dimensions of rates and reaction_coefficient_array, which will be nreactions-by-nspecies
        rates = np.swapaxes(rates, -2, -1)
        dc = rates @ self.reaction_coefficient_array

        # TODO: there will probably be an extra singleton dimension in here that we need to remove.
        return dc

    def do_timestep(self, dt, concentrations: np.ndarray, ambient_conditions: AmbientConditions, ambient_updated: bool) -> np.ndarray:
        # TODO: ODE solver
        # We should probably accept the solver as an init parameter so we can use scipy by default and switch to a custom one
        # (or a different scipy one) if needed.
        pass
