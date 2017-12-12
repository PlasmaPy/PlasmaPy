# PlasmaPy Development Roadmap

This document summarizes the overall development plans for PlasmaPy.
This roadmap is very fluid and is subject to change due to community
priorities.  For the most up-to-date information, see PlasmaPy's
[issues](https://github.com/PlasmaPy/plasmapy/issues), [pull
requests](https://github.com/PlasmaPy/plasmapy/pulls), and
[projects](https://github.com/PlasmaPy/PlasmaPy/projects).

## Features to be included in Version 0.1 release

This is to be a preview/prototype release.

- Plasma parameter calculations
- Transport coefficient functions
   - Braginskii theory
- Base PlasmaBlob class

## Features to be included in Version 0.2 release

- MHD simulation capability (finite difference & spectral)
- PIC simulation capability (only 1D for v0.1?)
- Simulation visualization/analysis tools
  - Basic turbulence analysis functionality from
    [TurbPlasma](https://github.com/tulasinandan/TurbPlasma) (to be
    speeded up via Numba/Cython)
    
## Features to be included in Version 0.3 release or at an indefinite time in future
- Magnetic topology analysis tools (data cube)
- Grad-Shafranov solver
- Dispersion solver
- Magnetic topology analysis tools (cylindrical geometry)
- MHD simulation capability (finite element)
- PIC simulation capability (2D)
- Spacecraft data analysis
- Experimental analysis tools
- Analytical plasma physics tools (with SymPy)
- Capability to read in CDF files (separate package?)

## Future events

- Code development meeting
- Software Carpentry workshops at plasma physics conferences
- Python in Plasma Physics conference?
