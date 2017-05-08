from setuptools import setup
import plasmapy

NAME = 'plasmapy'

"""The versioning scheme for PlasmaPy will be consistent with the
canonical format described extensively in PEP 440 [1] in order for
public versions to be compatible with setuptools and PyPI.  Releases
will use a [major].[minor].[micro] format with our first micro release
being '0.0.1' and our first minor release being '0.1'.  Development,
alpha, and beta versions will be appended by the pre-release tags
'.dev', 'a', and 'b', respectively, followed by an integer.

Here is a sample progression from older to newer versions: 
 
0.0.1.dev1  # early development version prior to 0.0.1
0.0.1.dev2
0.0.1a1     # alpha version prior to 0.0.1
0.0.1a2
0.0.1b1     # beta version prior to 0.0.1
0.0.1b2
0.0.1rc1    # release candidate
0.0.1       # micro release
0.1.dev1    # development version for first minor release
0.1a1       # alpha release prior to 0.1
0.1b1       # beta release prior to 0.1
0.1         # first minor release

The decision on when to change between development, alpha, and beta
versions will be made by the Code Development Committee with general
community agreement.  The trailing integer following pre-release tags
should be incremented either every one to two months, or more
frequently if appreciable changes to the code have accumulated.

[1] https://www.python.org/dev/peps/pep-0440/

"""

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
