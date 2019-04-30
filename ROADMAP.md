# PlasmaPy Development Roadmap

This document summarizes the long-term development plans for PlasmaPy.
This roadmap is very fluid and will continue to be updated and changed
throughout the development process based on community input.  

If you have feature requests, we highly encourage you to [raise an
issue](https://github.com/PlasmaPy/PlasmaPy/issues/new) in [PlasmaPy's
GitHub repository](https://github.com/PlasmaPy/PlasmaPy).  More
in-depth development plans are contained within the repository for
[PlasmaPy Enhancement Proposals
(PLEPs)](https://github.com/PlasmaPy/PlasmaPy-PLEPs) for some topics.
Additional resources (including grant proposals) may found on the
[PlasmaPy Community on
Zenodo](https://zenodo.org/communities/plasmapy).

## Objective

The goal of PlasmaPy is to foster a fully open source software
ecosystem for plasma research and education.  Functionality needed by
most plasma physicists will go in the PlasmaPy core package, whereas
more specialized functionality will go into PlasmaPy's affiliated
packages.

## Code development priorities

### High level data structures

One of the most fundamental design decisions for PlasmaPy is how to
represent different plasmas.  This decision is challenging given that
plasma data comes from many different sources, in many different
formats, and with many different conventions for storing and
representing metadata.  Creating a single base class will become
monolithic very quickly, whereas creating numerous different classes
has the potential problems of requiring users to remember many
different classes with diverging application programming interfaces
(APIs).  PlasmaPy is taking a hybrid approach (described in [PLEP
6](http://doi.org/10.5281/zenodo.1460977)) that involves creating an
[abstract base class
(ABC)](https://docs.python.org/3/library/abc.html) that contains
abstract methods that have to be defined for each of the subclasses
for different plasma representations (thus ensuring a more consistent
API).  Users can then invoke the `Plasma` class factory to select and
instantiate the appropriate subclass based on inputs provided by the
user.  The framework for base data structures was implemented by
Google Summer of Code student Ritiek Malhotra in 2018, along with the
[OpenPMD](https://github.com/openPMD/openPMD-standard) subclass.  The
next major development tasks are to create subclasses for different
plasma representations and to expand/refine the API defined by the
ABC.  

### Experimental support


### Plasma simulation capabilities


### Plasma parameters

- Dimensionless numbers functionality

- Code and documentation infrastructure
- Plasma parameter calculations
- Transport coefficient functions
   - Braginskii theory
- Base PlasmaBlob class


### Features to be included in Version 0.2 release

- MHD simulation capability (finite difference & spectral)
- PIC simulation capability (only 1D for v0.1?)
- Simulation visualization/analysis tools
  - Basic turbulence analysis functionality from
    [TurbPlasma](https://github.com/tulasinandan/TurbPlasma) (to be
    speeded up via Numba/Cython)
    
### Features to be included in Version 0.3 release or at an indefinite time in future
- Magnetic topology analysis tools (data cube)
- Grad-Shafranov solver
- Dispersion solver
- Magnetic topology analysis tools (cylindrical geometry)
- MHD simulation capability (finite element)
- PIC simulation capability (2D)
- Spacecraft data analysis
- Experimental analysis tools
- Analytical plasma physics tools (with SymPy)


### Future events

- Code development meeting
- Software Carpentry workshops at plasma physics conferences
- Python in Plasma Physics conference?
