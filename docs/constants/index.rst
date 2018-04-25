.. _plasmapy-constants:

********************************
Constants (`plasmapy.constants`)
********************************

The `~plasmapy.constants` module contains many physical constants that
are commonly used in plasma science.  With the exception of
`~plasmapy.constants.pi`, these constants are imported directly from
`~astropy.constants`.  We therefore recommend reviewing
`Astropy's constants subpackage documentation
<http://docs.astropy.org/en/stable/constants/index.html>`_.

.. warning::
    The values of some constants may be refined as better measurements
    become available.

    The values of some constants are not known precisely, and may change
    depending on the


    Units such as ``u.M_sun`` will use the current version of the
    corresponding constant. When using prior versions of the constants,
    quantities should be constructed with constants instead of units.
