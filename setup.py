from setuptools import setup
import plasmapy

NAME = 'plasmapy'

# The version scheme is being discussed at: 
#    https://github.com/PlasmaPy/PlasmaPy/issues/25

VERSION = '0.0.1.dev1'

setup(name=NAME,
      version=VERSION,
      description="Python package for plasma physics",
      requires=['numpy', 'scipy', 'astropy'],
      provides=[NAME],
      author="The PlasmaPy Community",
      author_email="namurphy@cfa.harvard.edu",  # until we get an email address
      license="BSD",
      url="https://github.com/PlasmaPy/PlasmaPy",  # until we make a webpage
      long_description=plasmapy.__doc__,
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
      packages=["plasmapy"],
      zip_safe=False,
      use_2to3=False,
      python_requires='>=3.6',
      )
