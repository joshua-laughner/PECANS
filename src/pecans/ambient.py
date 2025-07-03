from dataclasses import dataclass


@dataclass
class AmbientConditions:
    pressure: float
    temperature: float
    solar_zenith_angle: float