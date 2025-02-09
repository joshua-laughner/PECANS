from .grid_def import GridDef

class Pecans:
    def __init__(
        self,
        grid_definition: GridDef,
        transport_scheme,
        chemistry_scheme,
        initial_conditions,
        boundary_conditions,
        ambient_conditions, # temperature, pressure, sunlight
        winds,
        turbulence,
        emissions,
        output,
    ):
        self.grid_definition = grid_definition
        self.transport_scheme = transport_scheme
        self.chemistry_scheme = chemistry_scheme
        self.initial_conditions = initial_conditions
        self.boundary_conditions = boundary_conditions
        self.ambient_conditions = ambient_conditions
        self.winds = winds
        self.turbulence = turbulence
        self.emissions = emissions
        self.output = output