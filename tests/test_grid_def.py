import numpy as np
import pytest
from pecans.grid_def import RegularRectGrid, VariableRectGrid

# ------------------ #
# Regular grid tests #
# ------------------ #

@pytest.fixture
def reg_grid_1d():
    return RegularRectGrid(nx=10, ny=0, nz=0, dx=1000.0, dy=2000.0, dz=500.0)

@pytest.fixture
def reg_grid_2d():
    return RegularRectGrid(nx=10, ny=5, nz=0, dx=1000.0, dy=2000.0, dz=500.0)

@pytest.fixture
def reg_grid_3d():
    return RegularRectGrid(nx=10, ny=5, nz=3, dx=1000.0, dy=2000.0, dz=500.0)


class Test1DRegularRectGrid:
    def test_dims(self, reg_grid_1d):
        assert reg_grid_1d.grid_dims() == [10]

    def test_lengths(self, reg_grid_1d):
        lengths = reg_grid_1d.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.full(10, 1000.0))
        assert np.allclose(lengths[1], np.array([2000.0]))
        assert np.allclose(lengths[2], np.array([500.0]))

    def test_areas(self, reg_grid_1d):
        areas = reg_grid_1d.cell_areas()
        assert areas.shape == (10,)
        assert np.allclose(areas, np.full(10, 1000.0*2000.0))


class Test2DRegularRectGrid:
    def test_dims(self, reg_grid_2d):
        assert reg_grid_2d.grid_dims() == [10, 5]

    def test_lengths(self, reg_grid_2d):
        lengths = reg_grid_2d.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.full(10, 1000.0))
        assert np.allclose(lengths[1], np.full(5, 2000.0))
        assert np.allclose(lengths[2], np.array([500.0]))

    def test_areas(self, reg_grid_2d):
        areas = reg_grid_2d.cell_areas()
        assert areas.shape == (10,5)
        assert np.allclose(areas, np.full((10, 5), 1000.0*2000.0))


class Test3DRegularRectGrid:
    def test_dims(self, reg_grid_3d):
        assert reg_grid_3d.grid_dims() == [10, 5, 3]

    def test_lengths(self, reg_grid_3d):
        lengths = reg_grid_3d.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.full(10, 1000.0))
        assert np.allclose(lengths[1], np.full(5, 2000.0))
        assert np.allclose(lengths[2], np.full(3, 500.0))

    def test_areas(self, reg_grid_3d):
        areas = reg_grid_3d.cell_areas()
        assert areas.shape == (10,5)
        assert np.allclose(areas, np.full((10, 5), 1000.0*2000.0))


class TestAnyDRegularRectGrid:
    def test_invalid_shape(self):
        d = {'dx': 1.0, 'dy': 1.0, 'dz': 1.0}
        with pytest.raises(ValueError):
            # nx cannot be 0
            RegularRectGrid(nx=0, ny=0, nz=0, **d)
        with pytest.raises(ValueError):
            # no n can be negative
            RegularRectGrid(nx=1, ny=-1, nz=-1, **d)
        with pytest.raises(ValueError):
            # no n can be negative
            RegularRectGrid(nx=1, ny=1, nz=-1, **d)
        with pytest.raises(ValueError):
            # ny cannot be 0 if nz is not
            RegularRectGrid(nx=1, ny=0, nz=5, **d)

    def test_invalid_lengths(self):
        n = {'nx': 10, 'ny': 0, 'nz': 0}
        with pytest.raises(ValueError):
            # no length can be 0
            RegularRectGrid(**n, dx=1000.0, dy=0.0, dz=0.0)

        with pytest.raises(ValueError):
            # no length can be negative (checking dx)
            RegularRectGrid(**n, dx=-1000.0, dy=100.0, dz=200.0)
        with pytest.raises(ValueError):
            # no length can be negative (checking dy)
            RegularRectGrid(**n, dx=1000.0, dy=-100.0, dz=200.0)
        with pytest.raises(ValueError):
            # no length can be negative (checking dz)
            RegularRectGrid(**n, dx=1000.0, dy=100.0, dz=-200.0)

        with pytest.raises(ValueError):
            # no length can be NaN (checking dx)
            RegularRectGrid(**n, dx=np.nan, dy=100.0, dz=200.0)
        with pytest.raises(ValueError):
            # no length can be NaN (checking dy)
            RegularRectGrid(**n, dx=100.0, dy=np.nan, dz=200.0)
        with pytest.raises(ValueError):
            # no length can be NaN (checking dz)
            RegularRectGrid(**n, dx=200.0, dy=100.0, dz=np.nan)

        with pytest.raises(ValueError):
            # no length can be infinite (checking dx)
            RegularRectGrid(**n, dx=np.inf, dy=100.0, dz=200.0)
        with pytest.raises(ValueError):
            # no length can be NaN (checking dy)
            RegularRectGrid(**n, dx=100.0, dy=np.inf, dz=200.0)
        with pytest.raises(ValueError):
            # no length can be NaN (checking dz)
            RegularRectGrid(**n, dx=200.0, dy=100.0, dz=np.inf)


# ------------------------- #
# Variable width grid tests #
# ------------------------- #

@pytest.fixture
def var_grid_1d():
    return VariableRectGrid(
        dx=np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]),
        dy=np.array([5000.0]),
        dz=np.array([750.0])
    )

@pytest.fixture
def var_grid_2d():
    return VariableRectGrid(
        dx=np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]),
        dy=np.array([5000.0, 10000.0, 20000.0]),
        dz=np.array([750.0])
    )

@pytest.fixture
def var_grid_2d_one_y():
    return VariableRectGrid(
        dx=np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]),
        dy=np.array([5000.0]),
        dz=np.array([750.0]),
        n_dimensions=2,
    )

@pytest.fixture
def var_grid_3d():
    return VariableRectGrid(
        dx=np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]),
        dy=np.array([5000.0, 10000.0, 20000.0]),
        dz=np.array([750.0, 900.0])
    )

@pytest.fixture
def var_grid_3d_one_y():
    return VariableRectGrid(
        dx=np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]),
        dy=np.array([5000.0]),
        dz=np.array([750.0]),
        n_dimensions=3,
    )


class Test1DVariableRectGrid:
    def test_dims(self, var_grid_1d):
        assert var_grid_1d.grid_dims() == [5]

    def test_lengths(self, var_grid_1d):
        lengths = var_grid_1d.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]))
        assert np.allclose(lengths[1], np.array([5000.0]))
        assert np.allclose(lengths[2], np.array([750.0]))

    def test_areas(self, var_grid_1d):
        areas = var_grid_1d.cell_areas()
        assert areas.shape == (5,)
        assert np.allclose(areas, np.array([2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0]))


class Test2DVariableRectGrid:
    def test_dims(self, var_grid_2d):
        assert var_grid_2d.grid_dims() == [5, 3]

    def test_lengths(self, var_grid_2d):
        lengths = var_grid_2d.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]))
        assert np.allclose(lengths[1], np.array([5000.0, 10000.0, 20000.0]))
        assert np.allclose(lengths[2], np.array([750.0]))

    def test_areas(self, var_grid_2d):
        areas = var_grid_2d.cell_areas()
        assert areas.shape == (5,3)
        expected = np.array([
            [2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0],
            [2000.0*10000.0, 1000.0*10000.0, 500.0*10000.0, 250.0*10000.0, 125.0*10000.0],
            [2000.0*20000.0, 1000.0*20000.0, 500.0*20000.0, 250.0*20000.0, 125.0*20000.0],
        ]).T
        assert np.allclose(areas, expected)


class Test2DOneYVariableRectGrid:
    def test_dims(self, var_grid_2d_one_y):
        assert var_grid_2d_one_y.grid_dims() == [5,1]

    def test_lengths(self, var_grid_2d_one_y):
        lengths = var_grid_2d_one_y.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]))
        assert np.allclose(lengths[1], np.array([5000.0]))
        assert np.allclose(lengths[2], np.array([750.0]))

    def test_areas(self, var_grid_2d_one_y):
        areas = var_grid_2d_one_y.cell_areas()
        assert areas.shape == (5,1)
        assert np.allclose(
            areas,
            np.array([
                [2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0],
                [2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0],
            ]).T
        )


class Test3DVariableRectGrid:
    def test_dims(self, var_grid_3d):
        assert var_grid_3d.grid_dims() == [5, 3, 2]

    def test_lengths(self, var_grid_3d):
        lengths = var_grid_3d.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]))
        assert np.allclose(lengths[1], np.array([5000.0, 10000.0, 20000.0]))
        assert np.allclose(lengths[2], np.array([750.0, 900.0]))

    def test_areas(self, var_grid_3d):
        areas = var_grid_3d.cell_areas()
        assert areas.shape == (5,3)
        expected = np.array([
            [2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0],
            [2000.0*10000.0, 1000.0*10000.0, 500.0*10000.0, 250.0*10000.0, 125.0*10000.0],
            [2000.0*20000.0, 1000.0*20000.0, 500.0*20000.0, 250.0*20000.0, 125.0*20000.0],
        ]).T
        assert np.allclose(areas, expected)


class Test3DOneYVariableRectGrid:
    def test_dims(self, var_grid_3d_one_y):
        assert var_grid_3d_one_y.grid_dims() == [5,1,1]

    def test_lengths(self, var_grid_3d_one_y):
        lengths = var_grid_3d_one_y.cell_edge_lengths()
        assert len(lengths) == 3
        assert np.allclose(lengths[0], np.array([2000.0, 1000.0, 500.0, 250.0, 125.0]))
        assert np.allclose(lengths[1], np.array([5000.0]))
        assert np.allclose(lengths[2], np.array([750.0]))

    def test_areas(self, var_grid_3d_one_y):
        areas = var_grid_3d_one_y.cell_areas()
        assert areas.shape == (5,1)
        assert np.allclose(
            areas,
            np.array([
                [2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0],
                [2000.0*5000.0, 1000.0*5000.0, 500.0*5000.0, 250.0*5000.0, 125.0*5000.0],
            ]).T
        )


class TestAnyDVariableRectGrid:
    def test_non1d_lengths(self):
        # check that it catches if any axis is given as a 2D array
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones([5,3]), dy=np.ones(2), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones([5,3]), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones(2), dz=np.ones([5,3]))

    def test_empty_lengths(self):
        # Check if any axis is given as an empty array
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.array([]), dy=np.ones(2), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.array([]), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones(2), dz=np.array([]))

    def test_invalid_lengths(self):
        # Check if any axis has a negative or 0 length
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.array([-1.0]), dy=np.ones(2), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.array([0.0]), dy=np.ones(2), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.array([-1.0]), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.array([0.0]), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones(2), dz=np.array([-1.0]))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones(2), dz=np.array([0.0]))

        # Check if any axis has a NaN or infinite length
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.array([np.nan]), dy=np.ones(2), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.array([np.inf]), dy=np.ones(2), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.array([np.nan]), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.array([np.inf]), dz=np.ones(2))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones(2), dz=np.array([np.nan]))
        with pytest.raises(ValueError):
            VariableRectGrid(dx=np.ones(2), dy=np.ones(2), dz=np.array([np.inf]))
        