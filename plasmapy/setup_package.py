# script for astropy helpers to use to customize setup.
from setuptools.extension import Extension

# list of cython extensions in package
# cython_extensions = [Extension("*", ["*.pyx"])]

# list of c extensions in package
# c_extensions = [Extension("*", ["*.c"])]

# getting C and Cython extensions
cython_testFile = 'plasmapy/physics/parameters_cython.pyx'
c_testFile = 'plasmapy/physics/parameters_cython.c'

# names for compiling the .c files to
compileName = 'parameters_cython'
#compileNamePath = 'plasmapy/physics/parameters_cython'
compileNamePath = 'plasmapy.physics.parameters_cython'
compileName_cython = 'plasmapy/physics/parameters_cython_pyx'
compileName_c = 'plasmapy/physics/parameters_cython_c'


def get_extensions():
    extensions = []
    # getting cython files
    
    # getting .c files

    # expicitly adding cython and c file
    extensions.append(Extension(name=compileNamePath,
                                sources=[cython_testFile],
                                include_dirs=['numpy']))
#    extensions.append(Extension(name=compileName,
#                                sources=[c_testFile],
#                                include_dirs=['numpy']))
    return extensions
