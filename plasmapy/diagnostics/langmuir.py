"""
Defines the Langmuir analysis module as part of the diagnostics package.
"""

import numpy as np
import astropy.units as u


class Probe:
    r"""
    Class containing various Langmuir probe properties.

    For the moment cylindrical probes are assumed as spherical ones are
    generally impractical to make.

    Attributes
    ----------
    area : ~astropy.units.Quantity or ~numpy.ndarray
        Area or array of areas of the probe(s) in units convertible to m^2.

    config : str, optional
        Electrode configuration. Possible options are `'single'`
        (default), `'double'`, `'double-asym'`, `'triple'`, `'tetra'`,
        `'penta'`. The electrode configuration must correspond to the
        number of areas provided.

    """
    def __init__(self, area=u.m**2, config='single'):
        # Initialize probe properties
        self.area = area
        self.config = config


def swept_probe_analysis(potential_sweep,
                         measured_current,
                         probe,
                         ion='p+',
                         method='maxwellian',
                         smoothing=None,
                         polyorder=3):
    r"""
    Perform a Langmuir analysis of a given probe V-I profile in order to
    obtain various plasma parameters.

    Parameters
    ----------
    potential_sweep : ~astropy.units.Quantity
        The swept electric potential applied to the probe in units
        convertible to V.

    measured_current : ~astropy.units.Quantity
        The corresponding measured current in units convertible to A.

    probe : ~plasmapy.diagnostics.langmuir.Probe
        An instance of the Probe class containing probe properties.

    ion : str, optional
        Representation of the ion species. Defaults to Hydrogen.

    method : str, optional
        The preferred analysis method. Options are 'maxwellian' (default) and
        'OML', which applies Orbital Motion Limit theory.

    smoothing : str, optional
        The smoothing method applied to the measured current in order to
        obtain clean derivatives. Options are None (default) and 'savgol',
        which uses scipy.signal.savgol_filter. When `smoothing` is 'savgol',
        the polynomial order is given by polyorder.

    polyorder : int
        The order of the polynomial used to smooth the measured current and
        derivatives if `smoothing` is 'savgol'. Default is 3.

    Returns
    -------
    T_e : ~astropy.units.Quantity
        Best estimate of the electron temperature in units of eV.

    n_e : ~astropy.units.Quantity
        Best estimate of the electron density in units of m^-3.

    n_i : ~astropy.units.Quantity
        Best estimate of the ion density in units of m^-3. Equal to n_e
        unless the OML method is used.

    V_F : ~astropy.units.Quantity
        Best estimate of the floating potential in units of V.

    V_P : ~astropy.units.Quantity
        Best estimate of the plasma potential in units of V.

    Notes
    -----
    Due to the use of different methods to obtain `n_e` and `n_i`, the ratio
    between these does not necessarily give an accurate average ionization.
    """
    raise NotImplementedError


def obtain_EEDF(potential_sweep, measured_current, probe):
    r"""
    Obtain the Electron Energy Distribution Function of the plasma
    from the second derivative of the measured current response using
    Druyvestein's equations.

    Parameters
    ----------
    potential_sweep : ~astropy.units.Quantity
        The swept electric potential applied to the probe in units
        convertible to V.

    measured_current : ~astropy.units.Quantity
        The corresponding measured current in units convertible to A.

    probe : ~plasmapy.diagnostics.langmuir.Probe
        An instance of the Probe class containing probe properties.

    Returns
    -------
    EEDF : numpy.ndarray
        Array containing electron energies in units of eV and their
        corresponding probabilities in units of eV^-1 m^-3.
    """
    raise NotImplementedError


def obtain_EEPF(potential_sweep, measured_current, probe):
    r"""
    Obtain the Electron Energy Probability Function of the plasma
    from the second derivative of the measured current response using
    Druyvestein's equations.

    Parameters
    ----------
    potential_sweep : ~astropy.units.Quantity
        The swept electric potential applied to the probe in units
        convertible to V.

    measured_current : ~astropy.units.Quantity
        The corresponding measured current in units convertible to A.

    probe : ~plasmapy.diagnostics.langmuir.Probe
        An instance of the `~plasmapy.diagnostics.langmuir.Probe`
        class containing probe properties.

    Returns
    -------
    EEPF : numpy.ndarray
        Array of shape (2, N) containing electron energies in units of eV and
        their corresponding probabilities in units of eV^-3/2 m^-3.

    """
    raise NotImplementedError
