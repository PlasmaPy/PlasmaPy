from Cython.Build import cythonize
from setuptools import setup

setup(
    ext_modules=cythonize("plasmapy/formulary/speeds.pyx"),
    zip_safe=False,
)
