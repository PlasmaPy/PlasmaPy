# script for astropy helpers to use to customize setup.
from setuptools.extension import Extension

# list of cython extensions in package
# extensions = [Extension("*", ["*.pyx"])]

# getting C and Cython extensions
testFile = 'plasmapy/physics/parameters_cython.pyx'


def get_extensions():
    return [Extension(name='parameters_cython',
                      sources=[testFile],
                      include_dirs=['numpy'])]
