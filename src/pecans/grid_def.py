"""Classes to define the shape of the grid used by PECANS

The grid is the foundation of a box model; it determines how
finely the simulation region is discretized. While most other
parts of PECANS can be mixed and matched freely, some other
components will have to make assumptions about the grid shape.
The default components in PECANS will assume rectangular grid cells,
which have 2 neighbors per dimension. If you want to use a hexagonal
grid (for example), you should test each component to be sure it
works as expected with the new grid.
"""
from abc import ABC, abstractmethod
import numpy as np


class GridDef(ABC):
    """The base class for grid definitions.

    Any new grid definitions must implement the methods given here.
    """
    @abstractmethod
    def grid_dims(self) -> list[int]:
        """The list of lengths for each dimension of the grid.

        This must only return values for "active" dimensions. That is, for
        a 1D model with transport only along one dimension, this must return
        a list with length 1.
        """
        pass

    @abstractmethod
    def cell_edge_lengths(self) -> list[np.ndarray[tuple[int], np.dtype[np.float64]]]:
        """The lengths of each grid cell edge, in meters.

        The return value will be a list with one array per dimension;
        each array will have a number of elements equal to the model
        size in that dimension. Unlike ``grid_dims``, this will always
        return the maximum number of dimensions. For rectangular grids,
        this will be 3. As an example, a 2D model might return::

            [
              np.array([2000.0, 1000.0, 500.0, 1000.0, 200.0]),
              np.array([1000.0, 1000.0, 1000.0]),
              np.array([1000.0])
            ]

        The defines a 5-by-3 model domain where the boxes vary in length
        along the first dimension and have consistent length along the second.
        The third must still be defined, so that the boxes have a finite volume.
        """
        pass

    @abstractmethod
    def cell_areas(self) -> np.ndarray[tuple[int, ...], np.dtype[np.float64]]:
        """The areas of the bottom surface of the cells, in square meters.

        This must return an array with a shape matching one vertical layer of the model.
        For example, rectangular grids will return a 1D array for a 1D model and a 2D array
        for 2D or 3D models. This must report the area of the bottom surface of each grid
        cell, which is commonly needed for flux calculations.
        """
        pass


class RegularRectGrid(GridDef):
    """A grid where all grid cells are the same size and have two neighbors per dimension.

    This is the simplest grid. Each cell is a six-sided box and all cells are the same size.
    There are ``nx`` by ``ny`` by ``nz`` cells, with their sides' lengths given by ``dx``,
    ``dy``, and ``dz``. The lengths are given in meters.

    For a 1D model, set ``ny`` and ``nz`` to 0, and for a 2D model set ``nz`` to 0. ``nx``
    must always be >= 1. However, ``dx``, ``dy``, and ``dz`` must all be > 0, even in 1D or
    2D models, as the lengths may be needed for flux or other calculations.
    """
    def __init__(self, nx: int, ny: int, nz: int, dx: float, dy: float, dz: float):
        dims_check = {'nx': nx < 0, 'ny': ny < 0, 'nz': nz < 0}
        negative_dims = ', '.join(k for k, v in dims_check.items() if v)
        if negative_dims:
            raise ValueError(f'All dimensions must be >= 0, the following dim(s) were < 0: {negative_dims}')
        if nx < 1:
            raise ValueError('nx must be at least 1')
        if ny == 0 and nz > 0:
            raise ValueError('nz cannot be >0 if ny == 0')
        
        size_check = {'dx': dx <= 0.0, 'dy': dy <= 0.0, 'dz': dz <= 0.0}
        size_errors = ', '.join(k for k, v in size_check.items() if v)
        if size_errors:
            raise ValueError(f'All cell edge lengths must be > 0, the following lengths were <= 0: {size_errors}. '
                              'For a 1D or 2D model, all three lengths must still be nonzero; a reasonable default is to make the '
                              'unused dimension\'s lengths equal to one of the active dimensions (e.g., make dy = dx and dz = dx in a 1D model).')

        finite_check = {'dx': np.isfinite(dx), 'dy': np.isfinite(dy), 'dz': np.isfinite(dz)}
        finite_errors = ', '.join(k for k, v in finite_check.items() if not v)
        if finite_errors:
            raise ValueError(f'All cell edge lengths must be finite (not infinity or NaN), the following one(s) were not: {finite_errors}')
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def grid_dims(self) -> list[int]:
        if self.ny == 0:
            return [self.nx]
        elif self.nz == 0:
            return [self.nx, self.ny]
        else:
            return [self.nx, self.ny, self.nz]
    
    def cell_edge_lengths(self) -> list[np.ndarray[tuple[int], np.dtype[np.float64]]]:
        return [
            np.full(self.nx, self.dx),
            np.full(self.ny, self.dy),
            np.full(self.nz, self.dz),
        ]
    
    def cell_areas(self) -> np.ndarray[tuple[int, ...], np.dtype[np.float64]]:
        a = self.dx * self.dy
        if self.ny == 0:
            return np.full(self.nx, a)
        else:
            return np.full([self.nx, self.ny], a)
    

class VariableRectGrid(GridDef):
    """A grid where the grid cells' sizes may differ, but all cells still have two neighbors per dimension.

    This is similar to :class:`RegularRectGrid` in that it is a series of boxes, but in this case the boxes
    can have different lengths along a given dimension. These lengths are given by the arrays ``dx``, ``dy``,
    and ``dz``, which now must have a number of elements equal to the number of grid boxes along that dimension.
    For example, a 1D model with five boxes would use::

        dx = np.array([5000.0, 1000.0, 1000.0, 1000.0, 5000.0])
        dy = np.array([1000.0])
        dz = np.array([1000.0])

    This defines five boxes where the first and last are 5000 x 1000 x 1000 meters and the middle three are
    1000 x 1000 x 1000 meters. To define a 2D array with 3 boxes in the y direction, ``dy`` would have 3 elements,
    which could all be the same if you did not need to vary the box width. (:func:`np.full` will be useful in that case!)

    This grid will assume that it represents a 1D model if both ``dy`` and ``dz`` have length 1, and will assume it is
    2D if ``dz`` is one element long but ``dy`` has more than one element. To force a particular number of dimensions,
    use the ``n_dimensions`` argument. A value of ``1`` will ensure the model is 1D, a value of ``2`` will ensure the
    model is 2D. If the lengths of ``dx``, ``dy``, or ``dz`` are incompatible with that, a ``ValueError`` will be raised.
    """
    def __init__(
            self,
            dx: np.ndarray[tuple[int], np.dtype[np.float64]],
            dy: np.ndarray[tuple[int], np.dtype[np.float64]],
            dz: np.ndarray[tuple[int], np.dtype[np.float64]],
            n_dimensions: int | None = None
        ) -> None:
        arr_check = {'dx': np.ndim(dx), 'dy': np.ndim(dy), 'dz': np.ndim(dz)}
        bad_arrs = ', '.join(f'{k} is {v}D' for k, v in arr_check.items() if v != 1)
        if bad_arrs:
            raise ValueError(f'dx, dy, and dz must be 1D arrays, got {bad_arrs}')
        
        nx = np.size(dx)
        ny = np.size(dy)
        nz = np.size(dz)
        dims_check = {'nx': nx < 1, 'ny': ny < 1, 'nz': nz < 1}
        bad_dims = ', '.join(k for k, v in dims_check.items() if v)
        if bad_dims:
            raise ValueError(f'dx, dy, dz must each have at least one element, the following did not: {bad_dims}')
        
        size_check = {'dx': np.any(dx <= 0.0), 'dy': np.any(dy <= 0.0), 'dz': np.any(dz <= 0.0)}
        size_errors = ', '.join(k for k, v in size_check.items() if v)
        if size_errors:
            raise ValueError(f'All cell edge lengths must be > 0, the following lengths had elements <= 0: {size_errors}. '
                              'For a 1D or 2D model, all three lengths must still be nonzero; a reasonable default is to make the '
                              'unused dimension\'s lengths equal to one of the active dimensions (e.g., make dy = dx and dz = dx in a 1D model).')
        
        finite_check = {'dx': np.all(np.isfinite(dx)), 'dy': np.all(np.isfinite(dy)), 'dz': np.all(np.isfinite(dz))}
        finite_errors = ', '.join(k for k, v in finite_check.items() if not v)
        if finite_errors:
            raise ValueError(f'All elements of dx, dy, and dz must be finite (not infinities or NaNs), the following arrays had non-finite values: {finite_errors}')
        
        if n_dimensions is not None:
            if n_dimensions < 1:
                raise ValueError(f'n_dimensions must be >= 1, got {n_dimensions}')
            elif n_dimensions == 1:
                if ny > 1 or nz > 1:
                    raise ValueError('ny and/or nz cannot be >1 with n_dimensions = 1')
                ny, nz = 0, 0
            elif n_dimensions == 2:
                if nz > 1:
                    raise ValueError('nz cannot be >1 with n_dimensions = 2')
                nz = 0
        else:
            if ny == 1 and nz == 1:
                ny, nz = 0, 0
            elif nz == 1:
                nz = 0

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def grid_dims(self) -> list[int]:
        if self.ny == 0:
            return [self.nx]
        elif self.nz == 0:
            return [self.nx, self.ny]
        else:
            return [self.nx, self.ny, self.nz]
    
    def cell_edge_lengths(self) -> list[np.ndarray[tuple[int], np.dtype[np.float64]]]:
        return [self.dx, self.dy, self.dz]
    
    def cell_areas(self) -> np.ndarray[tuple[int, ...], np.dtype[np.float64]]:
        dx = self.dx.reshape(-1,1)
        dy = self.dy.reshape(1,-1)
        # Use squeeze to ensure that the 1D case does not return a 2D array
        # because of the broadcasting
        if self.ny == 0:
            return np.squeeze(dx * dy)
        else:
            return dx * dy
