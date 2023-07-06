"""
The `~plasmapy.dispersion.analytical` subpackage contains functionality
associated with analytical dispersion solutions.

.. attention::

   |expect-api-changes|
"""
__all__ = [
    "two_fluid",
    "stix",
    "AbstractMHDWave",
    "AlfvenWave",
    "FastMagnetosonicWave",
    "SlowMagnetosonicWave",
    "mhd_wave_class",
]

from plasmapy.dispersion.analytical.stix_ import stix
from plasmapy.dispersion.analytical.two_fluid_ import two_fluid
from plasmapy.dispersion.analytical.mhd_wave_class import (
    AbstractMHDWave,
    AlfvenWave,
    FastMagnetosonicWave,
    SlowMagnetosonicWave,
    mhd_waves,
)
