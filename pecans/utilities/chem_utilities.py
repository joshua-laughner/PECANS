from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np

from typing import Callable, Sequence


class ForcedParameter(ABC):
    @abstractmethod
    def current_value(self) -> np.ndarray:
        pass

    @abstractmethod
    def increment_time(self, dt):
        pass


@dataclass
class MechanismInterface:
    chem_solver: Callable[[float, dict, dict, dict], dict]
    species: Sequence[str]
    required_params: Sequence[str] = tuple()
    const_params: dict[str, np.ndarray] = field(default_factory=dict)
    forced_params: dict[str, ForcedParameter] = field(default_factory=dict)

    def call_solver(self, dt: float, species_in: dict):
        curr_forced_params = {k: fp.current_value() for k, fp in self.forced_params.items()}
        species_out = self.chem_solver(dt, self.const_params, curr_forced_params, species_in)
        for fp in self.forced_params.values():
            fp.increment_time(dt)
        return species_out
