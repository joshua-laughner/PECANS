from setuptools import setup
from Cython.Build import cythonize

setup(
    name="PECANS",
    #ext_modules=cythonize('pecans/chemderiv.pyx')  # accepts a glob pattern
    packages=['pecans']
)
