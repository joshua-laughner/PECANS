from abc import abstractmethod, ABC
import numpy as np


class AmbientConditions(ABC):
    @property
    @abstractmethod
    def pressure(self) -> np.ndarray:
        pass

    @property
    @abstractmethod
    def temperature(self) -> np.ndarray:
        pass


class MCMAmbient(AmbientConditions):
    def __init__(self, pressure: np.ndarray, temperature: np.ndarray, solar_zenith_angle: np.ndarray):
        self._pressure = pressure
        self._temperature = temperature
        self.solar_zenith_angle = solar_zenith_angle

    @property
    def pressure(self) -> np.ndarray:
        return self._pressure

    @property
    def temperature(self) -> np.ndarray:
        return self._temperature
