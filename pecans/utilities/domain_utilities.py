import numpy as np

def compute_coordinates_from_config(config):
    def coord_helper(n, d):
        # Since we start at the halfway point, ending at a full box will give the right number of boxes, since arange
        # excluded the stop value
        return np.arange(d / 2, n * d, d)

    dimensions = ['x', 'y', 'z']
    coords = []
    for dim in dimensions:
        d_dim = config.get('DOMAIN', 'd' + dim)
        n_dim = config.get('DOMAIN', 'n' + dim)
        if n_dim > 0:
            coords.append(coord_helper(n_dim, d_dim))
        else:
            coords.append(None)

    return tuple(coords)