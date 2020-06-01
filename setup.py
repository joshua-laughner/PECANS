from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
exts = [Extension(name='pecans.chemderiv', sources=['pecans/chemderiv.pyx'])]

setup(
    name="PECANS",
    ext_modules=cythonize(exts)  # accepts a glob pattern
    #packages=['pecans']
)
