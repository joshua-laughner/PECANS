

"""
Any init ideal methods need to return two values: the correct chem_solver function to call for each time step, which
must have the proper signiature (dt, TEMP, CAIR, species...), and a tuple of species names that should be initialized
for that model
"""

# TODO: make species definable elsewhere so that this could e.g. have multiple species with different emissions


def init_explicit_first_order_chem_solver(lifetime_seconds):
    lifetime_seconds = float(lifetime_seconds)

    def chem_solver(dt, TEMP, CAIR, **species):
        dt = float(dt)
        for specie, conc in species.items():
            # This ideal case assumes that all species are lost with the same first order rate constant, so
            # dC/dt = k*C, and k = 1/lifetime converted to seconds.
            species[specie] += -1 / lifetime_seconds * conc * dt

        return species

    return chem_solver, ('A',)
