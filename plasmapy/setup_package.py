# script for astropy helpers to use to customize setup.
from setuptools.extension import Extension

# list of cython extensions in package
# cython_extensions = [Extension("*", ["*.pyx"])]

# list of c extensions in package
# c_extensions = [Extension("*", ["*.c"])]

# getting C and Cython extensions
cython_testFile = 'plasmapy/physics/parameters_cython.pyx'
c_testFile = 'plasmapy/physics/parameters_cython.pyx'


def get_extensions():
    extensions = []
    # getting cython files
    
    # getting .c files

    # expicitly adding cython and c file
    extensions.append(Extension(name='parameters_cython.pyx',
                                sources=[cython_testFile],
                                include_dirs=['numpy']))
    extensions.append(Extension(name='parameters_cython.c',
                                sources=[c_testFile],
                                include_dirs=['numpy']))
    return extensions
