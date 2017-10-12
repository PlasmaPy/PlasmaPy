import ah_bootstrap
from setuptools import setup, find_packages
import builtins
builtins._ASTROPY_SETUP_ = True

#import astropy_helpers

from astropy_helpers.setup_helpers import (
    register_commands, get_package_info, get_debug_option)

# Package metadata
metadata = {}
with open('plasmapy/_metadata.py', 'r') as metadata_file:
    exec(metadata_file.read(), metadata)

# Requirements
with open('requirements/base.txt', 'r') as req_file:
    requirements = req_file.read().splitlines()

setup(name=metadata['name'],
      version=metadata['version'],
      description="Python package for plasma physics",
      requires=requirements,
      install_requires=requirements,
      provides=[metadata['name']],
      author=metadata['author'],
      author_email="namurphy@cfa.harvard.edu",  # until we get an email address
      license="BSD",
      url="https://github.com/PlasmaPy/PlasmaPy",  # until we make a webpage
      long_description=metadata['description'],
      keywords=['plasma', 'plasma physics', 'science'],
      classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Development Status :: 2 - Pre-Alpha',
        ],
      packages=find_packages(),
      zip_safe=False,
      use_2to3=False,
      python_requires='>=3.6',
      )
