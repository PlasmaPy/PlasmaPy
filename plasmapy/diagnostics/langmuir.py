"""
Defines the Langmuir analysis module as part of the diagnostics package.
"""

import plasmapy.constants as const
import astropy.units as u
import copy

from plasmapy import utils

from astropy.visualization import quantity_support
import numpy as np
import matplotlib.pyplot as plt

def fit_func_lin(x, x0, y0, c0):
    r"""Linear fitting function."""
    return y0 + c0 * (x - x0)

def fit_func_double_lin(x, x0, y0, c0, c1):
    r"""Piecewise linear fitting function. (x0, y0) denotes the location of
    the knee of the transition, with c0 and c1 being the left and right-hand
    slopes, respectively."""
    return np.piecewise(x, x < x0, [lambda x: y0 + c0 * (x - x0), 
                                    lambda x: y0 + c1 * (x - x0)])  
    
class Characteristic:
    r"""Class representing a I-V probe characteristic for convenient access.
    Supports units.

    Attributes
    ----------
    bias : Quantity, ndarray
        Array of applied probe biases in units convertible to V.

    current : Quantity, ndarray
        Array of applied probe currents in units convertible to A.
        
    """
    
    def __init__(self, bias, current):
        self.bias = bias
        self.current = current
        
        self.check_validity()

    def __getitem__(self, key):
        r"""Allows array indexing operations directly on the Characteristic
        object.
        
        """
        
        b = Characteristic(self.bias[key], self.current[key])
        return b

    def __sub__(self, other):
        b = copy.deepcopy(self)
        b.current -= other.current
        return b

    def __add__(self, other):
        b = copy.deepcopy(self)
        b.current += other.current
        return b
    
    def sort(self):
        r"""Sorts the characteristic into an ascending bias.
        
        """
        
        _sort = self.bias.argsort()
        self.current = self.current[_sort]
        self.bias = self.bias[_sort]
    
    def get_unique(self):
        r"""Removes duplicate bias values through averaging and returns this
        as a new object.
        
        """
        
        bias_unique = np.unique(self.bias)
        current_unique = []
        for bias in bias_unique:
            current_unique = np.append(current_unique, np.mean(self.current[self.bias == bias].to(u.A).value))
        current_unique *= u.A
        
        return Characteristic(bias_unique, current_unique)
    
    def check_validity(self):
        r"""Checks the unit and value validity of the characteristic.
        
        """
        
        utils._check_quantity(self.bias, 'bias', str(self), u.V,
                              can_be_negative=True, can_be_complex=False, 
                              can_be_inf=False, can_be_nan=True)
        
        utils._check_quantity(self.current, 'bias', str(self), u.A,
                              can_be_negative=True, can_be_complex=False, 
                              can_be_inf=False, can_be_nan=True)
        
    def get_padded_limit(self, padding, log=False):  
        r"""Returns the limits of the current range for plotting, taking into
        account padding. Matplotlib lacks this functionality.
        
        Parameters
        ----------
        padding : float
            The padding ratio as a float between 0.0 and 1.0.
            
        log : boolean, optional
            If True the calculation will be performed on a logarithmic scale. 
            Default is False.
        
        """
        
        ymax = np.max(self.current).to(u.A).value
        
        if(log):
            ymin = np.min(np.abs(self.current[self.current != 0])).to(u.A).value
            return [ymin * 10**(-padding * np.log10(ymax / ymin)), 
                    ymax * 10**(padding * np.log10(ymax / ymin))] * u.A
        else:  
            ymin = np.min(self.current).to(u.A).value
            return [ymin - padding * (ymax - ymin), ymax + padding * (ymax - ymin)] * u.A
        
    def plot(self):
        r"""Plots the characteristic in matplotlib.
        
        """
        
        with quantity_support(): 
            plt.figure()
            plt.scatter(self.bias, self.current, 
                        marker='.', color='k')
            plt.title("Probe characteristic")

@utils.check_quantity({'probe_area': {'units': u.m**2, 
                                      'can_be_negative': False, 
                                      'can_be_complex': False, 
                                      'can_be_inf': False, 
                                      'can_be_nan': False}})
def swept_probe_analysis(probe_characteristic, probe_area, gas, 
                         bimaxwellian=False, visualize=False,
                         plot_EEDF=False):
    r"""Attempts to perform a basic swept probe analysis based on the provided
    characteristic and probe data. Suitable for single cylindrical probes in
    low-pressure DC plasmas, since OML is applied.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The swept probe characteristic that is to be analyzed.

    probe_area : Quantity
        The area of the probe exposed to plasma in units convertible to m^2.

    gas : float
        The (mean) mass of the background gas in atomic mass units (amu).

    visualize : bool, Optional
        Can be used to plot the characteristic and the obtained parameters.
        Default is False.

    plot_EEDF : bool, Optional
        If True, the EEDF is computed and shown. Default is False.

    Returns
    -------
    T_e : Quantity
        Best estimate of the electron temperature in units of eV. Contains
        two values if bimaxwellian is True.

    n_e : Quantity
        Estimate of the electron density in units of m^-3. See the Notes on 
        plasma densities.

    n_i : Quantity
        Estimate of the ion density in units of m^-3. See the Notes on 
        plasma densities.

    n_i_OML : Quantity
        OML-theory estimate of the ion density in units of m^-3. See the Notes 
        on plasma densities.

    V_F : Quantity
        Estimate of the floating potential in units of V.

    V_P : Quantity
        Estimate of the plasma potential in units of V.

    I_e : Quantity
        Estimate of the electron saturation current in units of Am^-2.

    I_i : Quantity
        Estimate of the ion saturation current in units of Am^-2.

    hot_fraction : float
        Estimate of the total hot (energetic) electron fraction.

    Notes
    -----
    This function combines the separate probe analysis functions into a single
    analysis. Results are returned as a Dictionary. On plasma densities: in an 
    ideal quasi-neutral plasma all densities should be equal. However, in 
    practice this will not be the case. The electron density is the poorest
    estimate due to the hard to obtain knee in the electron current. The
    density provided by OML theory is likely the best estimate as it is not
    dependent on the obtained electron temperature, given that the conditions
    for OML theory hold.
    
    """
    
    # Check (unit) validity of the probe characteristic
    probe_characteristic.check_validity()
        
    if(visualize):
        with quantity_support():                
            fig, (ax1, ax2) = plt.subplots(2, 1)
            ax1.plot(probe_characteristic.bias, probe_characteristic.current, 
                    marker='.', color='k', linestyle='')
            ax1.set_title("Probe characteristic") 
            ax2.set_ylim(probe_characteristic.get_padded_limit(0.1))
            
            ax2.plot(probe_characteristic.bias, 
                     np.abs(probe_characteristic.current), \
                            marker='.', color='k', linestyle='')
            ax2.set_title("Logarithmic") 
            ax2.set_ylim(probe_characteristic.get_padded_limit(0.1, log=True))
            
    # Obtain the plasma and floating potentials
    V_P = get_plasma_potential(probe_characteristic)
    V_F = get_floating_potential(probe_characteristic)
    
    # Obtain the electron and ion saturation currents
    I_e = get_electron_saturation_current(probe_characteristic)
    I_i = get_ion_saturation_current(probe_characteristic)
    
    # The OML method is used to obtain an ion density without knowing the
    # electron temperature. This can then be used to obtain the ion current
    # and subsequently a better electron current fit.
    n_i_OML, fit = get_ion_density_OML(probe_characteristic, probe_area, gas,
                                       return_fit=True)
    
    ion_current = extrapolate_ion_current_OML(probe_characteristic, fit)
    
    # First electron temperature iteration
    exponential_section = extract_exponential_section(probe_characteristic, 
                                                      ion_current=ion_current)   
    T_e, hot_fraction = get_electron_temperature(exponential_section, 
                                      bimaxwellian=bimaxwellian, 
                                      return_hot_fraction=True)
    
    # Second electron temperature iteration, using an electron temperature-
    # adjusted exponential section
    exponential_section = extract_exponential_section(probe_characteristic, 
                                                      T_e=T_e,
                                                      ion_current=ion_current)       
    T_e, hot_fraction, fit = get_electron_temperature(exponential_section,
                                                      bimaxwellian=bimaxwellian,
                                                      return_fit=True, 
                                                      return_hot_fraction=True)
    
    # Extrapolate the fit of the exponential section to obtain the full 
    # electron current. This has no use in the analysis except for 
    # visualization.
    electron_current = extrapolate_electron_current(probe_characteristic,
                                                    fit,
                                                    bimaxwellian=bimaxwellian)
    
    # Using a good estimate of electron temperature, obtain the ion and
    # electron densities from the saturation currents.
    n_i = get_ion_density_LM(I_i, reduce_bimaxwellian_temperature(T_e, 
                             hot_fraction), probe_area, gas)
    n_e = get_electron_density_LM(I_e, reduce_bimaxwellian_temperature(T_e, 
                                  hot_fraction), probe_area)
    
    if(visualize):
        with quantity_support():  
            ax1.axvline(x=V_P.value, color='gray', linestyle='--')
            ax1.axhline(y=I_e.value, color='grey', linestyle='--')
            ax1.axvline(x=V_F.value, color='k', linestyle='--')
            ax1.axhline(y=I_i.value, color='r', linestyle='--')
            ax1.plot(ion_current.bias, ion_current.current, c='y')
            ax1.plot(electron_current.bias, electron_current.current, c='c')
            tot_current = ion_current + electron_current
            ax1.plot(tot_current.bias, tot_current.current, c='g')
            
            ax2.axvline(x=V_P.value, color='gray', linestyle='--')
            ax2.axhline(y=I_e.value, color='grey', linestyle='--')
            ax2.axvline(x=V_F.value, color='k', linestyle='--')
            ax2.axhline(y=I_i.value, color='r', linestyle='--')
            ax2.plot(ion_current.bias, np.abs(ion_current.current), c='y')
            ax2.plot(electron_current.bias, np.abs(electron_current.current), c='c')
            ax2.plot(tot_current.bias, np.abs(tot_current.current), c='g')
            ax2.set_yscale("log", nonposy='clip')
            
            fig.tight_layout()
         
    # Obtain and show the EEDF. This is only useful if the characteristic data
    # has been preprocessed to be sufficiently smooth and noiseless.
    if(plot_EEDF):
        get_EEDF(probe_characteristic, probe_area, n_i_OML, visualize=True)
    
    # Compile the results dictionary
    results = {'V_P': V_P, 'V_F': V_F, 'I_e': I_e, 'I_i': I_i, 'n_e': n_e, 
            'n_i': n_i, 'T_e': T_e, 'n_i_OML': n_i_OML}
    
    if(bimaxwellian):
        results['hot_fraction'] = hot_fraction
    
    return results
    
def get_plasma_potential(probe_characteristic, return_arg=False):
    r"""Implements the simplest but crudest method for obtaining an estimate of 
    the plasma potential from the probe characteristic.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.

    return_arg : bool, optional
        Controls whether or not the argument of the plasma potential within the
        characteristic array should be returned instead of the value of the 
        voltage. Default is False.

    Returns
    -------
    V_P : Quantity
        Estimate of the plasma potential in units convertible to V.

    Notes
    -----
    The method used in the function takes the maximum gradient of the probe 
    current as the 'knee' of the transition from exponential increase into the 
    electron the saturation region.
    
    """
    
    probe_characteristic.check_validity()
    
    # Sort the characteristic prior to differentiation
    probe_characteristic.sort()
    
    characteristic_unique = probe_characteristic.get_unique()
    
    # Crude method of acquiring the first derivative
    dIdV = np.diff(characteristic_unique.current.to(u.A).value) / \
           np.diff(characteristic_unique.bias.to(u.V).value)
         
    arg_V_P = np.argmax(dIdV)
        
    if(return_arg):
        return characteristic_unique.bias[arg_V_P], arg_V_P
    
    return characteristic_unique.bias[arg_V_P]

def get_floating_potential(probe_characteristic, return_arg=False):
    r"""Implements the simplest but crudest method for obtaining an estimate of 
    the floating potential from the probe characteristic.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.

    return_arg : bool, optional
        Controls whether or not the argument of the floating potential within 
        the characteristic array should be returned instead of the value of the 
        voltage. Default is False.

    Returns
    -------
    V_F : Quantity
        Estimate of the floating potential in units convertible to V.

    Notes
    -----
    The method used in this function takes the probe current closest to zero 
    Amperes as the floating potential.
    
    """
    
    probe_characteristic.check_validity()
    
    arg_V_F = np.argmin(np.abs(probe_characteristic.current))
    
    if(return_arg):
        return probe_characteristic.bias[arg_V_F], arg_V_F
    
    return probe_characteristic.bias[arg_V_F]

def get_electron_saturation_current(probe_characteristic):
    r"""Obtains an estimate of the electron saturation current corresponding 
    to the obtained plasma potential.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.

    Returns
    -------
    I_e : Quantity
        Estimate of the electron saturation current in units convertible to A.

    Notes
    -----
    The function `get_plasma_potential` is used to obtain an estimate of the
    plasma potential. The corresponding electron saturation current is 
    returned.
    
    """
    
    probe_characteristic.check_validity()
    
    _, arg_V_P = get_plasma_potential(probe_characteristic, return_arg=True)
    
    return probe_characteristic.current[arg_V_P]

def get_ion_saturation_current(probe_characteristic):
    r"""Implements the simplest but crudest method for obtaining an estimate of 
    the ion saturation current from the probe characteristic.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.

    Returns
    -------
    I_i : Quantity
        Estimate of the ion saturation current in units convertible to A.

    Notes
    -----
    The method implemented in this function takes the ion saturation current 
    as the smallest probe current in the characteristic. This assumes the bias 
    range in the ion region is sufficiently negative for the ion current to
    saturate.
    
    """
    
    probe_characteristic.check_validity()
    
    return np.min(probe_characteristic.current)

@utils.check_quantity({'probe_area': {'units': u.m**2, 
                                      'can_be_negative': False, 
                                      'can_be_complex': False, 
                                      'can_be_inf': False, 
                                      'can_be_nan': False}})
def get_ion_density_LM(ion_saturation_current, T_e, 
                       probe_area, gas):
    r"""Implements the Langmuir-Mottley (LM) method of obtaining the ion 
    density.

    Parameters
    ----------
    ion_saturation_current : Quantity
        The ion saturation current in units convertible to A.
        
    electron_temperature : Quantity
        The electron temperature in units convertible to eV.

    Returns
    -------
    n_i : Quantity
        Estimate of the ion density in units convertible to m^-3.

    Notes
    -----
    The method implemented in this function obtains the ion density from the
    ion saturation current density assuming that the ion current loss to the 
    probe is equal to the Bohm loss. The acoustic Bohm velocity is obtained
    from the electron temperature and the ion mass.
    
    """
    
    utils._check_quantity(T_e, 'T_e', 'extract_exponential_section', u.eV, 
                          can_be_negative=False)
    
    c_s = np.sqrt(T_e / (gas * const.m_n))
    
    n_i = np.abs(ion_saturation_current) / (const.e * probe_area * c_s)
    
    return n_i.to(u.m**-3)

@utils.check_quantity({'probe_area': {'units': u.m**2, 
                                      'can_be_negative': False, 
                                      'can_be_complex': False, 
                                      'can_be_inf': False, 
                                      'can_be_nan': False}})
def get_electron_density_LM(electron_saturation_current, T_e, 
                            probe_area):
    r"""Implements the Langmuir-Mottley (LM) method of obtaining the electron 
    density.

    Parameters
    ----------
    electron_saturation_current : Quantity
        The electron saturation current in units convertible to A.
        
    electron_temperature : Quantity
        The electron temperature in units convertible to eV.

    Returns
    -------
    n_e : Quantity
        Estimate of the electron density in units convertible to m^-3.

    Notes
    -----
    The method implemented in this function obtains the electron density from 
    the electron saturation current density, assuming a low plasma density. 
    Please note that the electron saturation current density is a hard
    parameter to acquire and it is usually better to measure the ion density,
    which should be identical to the electron density in quasineutral plasmas.
    
    """
    
    utils._check_quantity(T_e, 'T_e', 'extract_exponential_section', u.eV, 
                          can_be_negative=False)
     
    v_e = np.sqrt(T_e / (2 * np.pi * const.m_e))
    
    n_e = electron_saturation_current / (probe_area * const.e * v_e)
    
    return n_e.to(u.m**-3)

def extract_exponential_section(probe_characteristic, T_e=None, 
                                ion_current=None):
    r"""Extracts the section of exponential electron current growth from the
    probe characteristic.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.
        
    T_e : Quantity, optional
        The electron temperature can be supplied to improve the accuracy of the 
        bounds of the exponential region.

    Returns
    -------
    exponential_section : Characteristic
        The exponential electron current growth section.

    Notes
    -----
    This function extracts the region of exponential electron growth from the 
    probe characteristic under the assumption that this bias region is bounded
    by the floating and plasma potentials. Additionally, an improvement in 
    accuracy can be made when the electron temperature is supplied.
    
    """
    
    probe_characteristic.check_validity()
    
    V_F = get_floating_potential(probe_characteristic)
    
    V_P = get_plasma_potential(probe_characteristic)
    
    if(T_e != None):
        
        utils._check_quantity(T_e, 'T_e', 'extract_exponential_section', u.eV, 
                              can_be_negative=False)
        
        # If a bi-Maxwellian electron temperature is supplied grab the first 
        # (cold) temperature
        if(np.array(T_e).size > 1):
            T_e = np.min(T_e)
        
        _filter = (probe_characteristic.bias > V_F + 1.5 * T_e / const.e) & (
                probe_characteristic.bias < V_P - 0.2 * T_e / const.e)
    else:
        _filter = (probe_characteristic.bias > V_F) & (
                probe_characteristic.bias < V_P)
        
    exponential_section = probe_characteristic[_filter]
    
    if(ion_current != None):
        exponential_section = exponential_section - ion_current[_filter]
    
    return exponential_section

def extract_ion_section(probe_characteristic):
    r"""Extracts the section dominated by ion collection from the probe 
    characteristic.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.

    Returns
    -------
    ion_section : Characteristic
        The exponential electron current growth section.

    Notes
    -----
    This function extracts the region dominated by ion collection from the 
    probe characteristic under the assumption that this bias region is only
    bounded by the floating potential on the right hand side.
    
    """
    
    probe_characteristic.check_validity()
    
    V_F = get_floating_potential(probe_characteristic)
    
    ion_section = probe_characteristic[probe_characteristic.bias < V_F]
    
    return ion_section

def get_electron_temperature(exponential_section, visualize=False,
                             bimaxwellian=False, return_fit=False,
                             return_hot_fraction=False):
    r"""Extracts the section dominated by ion collection from the probe 
    characteristic.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.

    visualize : bool, optional
        If True a plot of the exponential fit is shown. Default is False.
        
    return_fit: bool, optional
        If True the parameters of the fit will be returned in addition to the 
        electron temperature. Default is False.
        
    return_hot_fraction: float, optional
        If True the function will return the hot fraction in case one is 
        obtained from a bi-Maxwellian fit.

    Returns
    -------
    T_e : Quantity
        The estimated electron temperature in eV.

    Notes
    -----
    The exponential section of the probe characteristic should be a straight 
    line if the plasma electrons are fully Maxwellian. The slope is then 
    inversely proportional to the electron temperature.
    
    """
    
    exponential_section.check_validity()
    
    # Remove values in the section with a current equal to or smaller than 
    # zero.
    exponential_section = exponential_section[
            exponential_section.current.to(u.A).value > 0]

    p0 = None
    
    # Instantiate the correct fitting equation
    if(bimaxwellian):
        x0 = np.min(exponential_section.bias) + 2/3 * \
                         (np.max(exponential_section.bias) - \
                          np.min(exponential_section.bias))
        
        p0 = [x0.to(u.V).value, 0.6, 0, 5]
        
        fit_func = fit_func_double_lin
    else:
        fit_func = fit_func_lin
            
    from scipy.optimize import curve_fit   
    
    # Perform the actual fit of the data
    fit, _ = curve_fit(fit_func, exponential_section.bias.to(u.V).value, 
                       np.log(exponential_section.current.to(u.A).value), 
                       p0=p0)
    
    hot_fraction = None
    
    # Obtain the plasma parameters from the fit
    if(not bimaxwellian):
        slope = fit[2]
        
        T_e = slope**-1 * u.eV
    else:
        x0, y0 = fit[0], fit[1]
        slope1, slope2 = [fit[2], fit[3]]
        
        # In order to obtain the energetic electron fraction the fits of the 
        # cold and hot populations are extrapolated to the plasma potential
        # (ie. the maximum bias of the exponential section). The logarithmic
        # difference between these currents equates to the density difference.
        k1 = fit_func_lin(np.max(exponential_section.bias.to(u.V).value),
                          *fit[[0,1,3]])
        k2 = fit_func_lin(np.max(exponential_section.bias.to(u.V).value),
                          *fit[0:3])
        
        # Compute the total hot (energetic) fraction
        hot_fraction = 1 / (1 + np.exp(k1 - k2))
        
        # If bi-Maxwellian, return main temperature first
        T_e = np.array([slope2**-1, slope1**-1]) * u.eV
        
    if(visualize):
        with quantity_support():   
            plt.figure()
            plt.scatter(exponential_section.bias.to(u.V), 
                     np.log(exponential_section.current.to(u.A).value), 
                     color='k', marker='.')
            if(bimaxwellian):
                plt.scatter(x0, y0, marker='o', c='g')
                plt.plot(exponential_section.bias.to(u.V), 
                         fit_func_lin(exponential_section.bias.to(u.V).value, 
                                      *fit[0:3]), c='g', linestyle='--')
            plt.plot(exponential_section.bias.to(u.V), 
                     fit_func(exponential_section.bias.to(u.V).value, *fit), 
                     c='g')
            plt.ylabel("Logarithmic current")
            plt.title("Exponential fit")
            plt.tight_layout()
        
    k = [T_e]
    
    if(return_hot_fraction):
        k.append(hot_fraction)
    
    if(return_fit):
        k.append(fit)
        
    return k

def extrapolate_electron_current(probe_characteristic, fit,
                                 bimaxwellian=False, visualize=False):
    r"""Extrapolates the electron current from the Maxwellian electron 
    temperature obtained in the exponential growth region.
    
    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.
        
    fit : ndarray
        Polynomial fit coefficients returned by the electron temperature fit.

    visualize : bool, optional
        If True a plot of the extracted electron current is shown. Default is 
        False.

    Returns
    -------
    electron_current : Characteristic
        The extrapolated electron current characteristic.

    Notes
    -----
    Assuming the electron population is fully Maxwellian the pure electron 
    current is extrapolated from the fit of the exponential region for the
    entire bias range.
    
    """
    
    probe_characteristic.check_validity()
        
    if(bimaxwellian):
        fit_func = fit_func_double_lin
    else:
        fit_func = fit_func_lin
        
    electron_current = np.exp(fit_func(probe_characteristic.bias.to(u.V).value, 
                                       *fit))*u.A
    
    electron_characteristic = Characteristic(probe_characteristic.bias, 
                                             electron_current)
    
    electron_characteristic.current[electron_characteristic.current > np.max(
            probe_characteristic.current)] = np.NaN
        
    if(visualize):
        with quantity_support():   
            plt.figure()
            plt.scatter(probe_characteristic.bias, 
                        probe_characteristic.current.to(u.mA), marker='.', 
                        c='k')
            plt.plot(electron_characteristic.bias, 
                     electron_characteristic.current.to(u.mA))
    
    return electron_characteristic

@utils.check_quantity({'T_e': {'units': u.eV, 
                               'can_be_negative': False, 
                               'can_be_complex': False, 
                               'can_be_inf': False, 
                               'can_be_nan': False}})
def reduce_bimaxwellian_temperature(T_e, hot_fraction):
    r"""Function to reduce a bi-Maxwellian (dual) temperature to a single 
    mean temperature for a given fraction.

    Parameters
    ----------
    T_e : Quantity, ndarray
        The bi-Maxwellian temperatures in eV. If a single temperature is 
        provided, this is returned.
        
    hot_fraction : float
        Fraction of hot to total population. If this parameter is None the 
        temperature is assumed to be singular Maxwellian.

    Returns
    -------
    T_e : Quantity
        The reduced (mean) temperature in units of eV.

    Notes
    -----
    This function aids methods that take a single electron temperature in 
    situations where the electron population is bi-Maxwellian. The reduced
    temperature is obtained as the weighted mean.
    
    """
    
    # Return the electron temperature itself if it is not bi-Maxwellian
    # in the first place.
    if(hot_fraction == None or not np.array(T_e).size > 1):
        return T_e
    
    return T_e[0] * (1 - hot_fraction) + T_e[1] * hot_fraction

@utils.check_quantity({'probe_area': {'units': u.m**2, 
                                      'can_be_negative': False, 
                                      'can_be_complex': False, 
                                      'can_be_inf': False, 
                                      'can_be_nan': False}})
def get_ion_density_OML(probe_characteristic, probe_area, gas, 
                        visualize=False, return_fit=False):
    r"""Implements the Orbital Motion Limit (OML) method of obtaining an
    estimate of the ion density.

    Parameters
    ----------
    probe_characteristic : Characteristic
        The swept probe characteristic that is to be analyzed.

    probe_area : Quantity
        The area of the probe exposed to plasma in units convertible to m^2.

    gas : float
        The (mean) mass of the background gas in atomic mass units (amu).

    visualize : bool, optional
        If True a plot of the OML fit is shown. Default is False.
        
    return_fit: bool, optional
        If True the parameters of the fit will be returned in addition to the 
        ion density. Default is False.

    Returns
    -------
    n_i_OML : Quantity
        Estimated ion density in m^-3.

    Notes
    -----
    The method implemented in this function holds for cylindrical probes in a 
    cold ion plasma, ie. T_i = 0 eV. With OML theory an expression is found
    for the ion current as function of probe bias independent of the electron 
    temperature:
    
    .. math::
        I = A_p n e \frac{\sqrt{2}}{\pi} \left( \frac{|eV_b|}{m_i} 
        \right)^{1/2} 
        
    """
    
    probe_characteristic.check_validity()
    
    ion_section = extract_ion_section(probe_characteristic)
    
    fit = np.polyfit(ion_section.bias.to(u.V).value,
                   ion_section.current.to(u.mA).value**2,
                   1)
    
    poly = np.poly1d(fit)
    
    slope = fit[0]
    
    n_i_OML = np.sqrt(-slope*u.mA**2/u.V * np.pi**2 * (gas * const.m_n) / 
                      (probe_area**2 * const.e**3 * 2))
                            
    if(visualize):
        with quantity_support():   
            plt.figure()
            plt.scatter(ion_section.bias.to(u.V), ion_section.current.to(
                    u.mA)**2, 
                     color='k', marker='.')
            plt.plot(ion_section.bias.to(u.V), 
                     poly(ion_section.bias.to(u.V).value), 
                     c='g')
            plt.title("OML fit")
            plt.tight_layout()
                
    if(return_fit):
        return n_i_OML.to(u.m**-3), fit
    
    return n_i_OML.to(u.m**-3)

def extrapolate_ion_current_OML(probe_characteristic, fit, 
                            visualize=False):
    r"""Extrapolates the ion current from the ion density obtained with the
    OML method.
        
    Parameters
    ----------
    probe_characteristic : Characteristic
        The probe characteristic that is being analyzed.
        
    fit : ndarray
        Fit coefficients returned by the OML method.

    visualize : bool, optional
        If True a plot of the extracted electron current is shown. Default is 
        False.

    Returns
    -------
    ion_section : Characteristic
        The exponential electron current growth section.

    Notes
    -----
    The exponential section of the probe characteristic should be a straight 
    line if the plasma electrons are fully Maxwellian. The slope is then 
    inversely proportional to the electron temperature.
    
    """
    
    probe_characteristic.check_validity()
    
    slope = fit[0] * u.mA**2 / u.V
    offset = fit[1] * u.mA**2
    
    ion_current = -np.sqrt(np.clip(slope * probe_characteristic.bias + offset, 
                                   0.0, None))
    
    ion_characteristic = Characteristic(probe_characteristic.bias, ion_current)
    
    if(visualize):
        with quantity_support():   
            plt.figure()
            plt.scatter(probe_characteristic.bias, 
                        probe_characteristic.current.to(u.mA), marker='.', 
                        c='k')
            plt.plot(probe_characteristic.bias, 
                     ion_characteristic.current.to(u.mA), c='y')
    
    return ion_characteristic
    
@utils.check_quantity({'probe_area': {'units': u.m**2, 
                                      'can_be_negative': False, 
                                      'can_be_complex': False, 
                                      'can_be_inf': False, 
                                      'can_be_nan': False}})
def get_EEDF(probe_characteristic, probe_area, n_i, visualize=False):
    r"""Implements the Druyvesteyn method of obtaining the normalized 
    Electron Energy Distribution Function (EEDF).

    Parameters
    ----------
    probe_characteristic : Characteristic
        The swept probe characteristic that is to be analyzed.

    probe_area : Quantity
        The area of the probe exposed to plasma in units convertible to m^2.

    n_i : Quantity
        The ion density in units covertible to m^-3.

    visualize : bool, optional
        If True a plot of the extracted electron current is shown. Default is 
        False.

    Returns
    -------
    energy : Quantity, ndarray
        Array of potentials in V.
        
    probability : float, ndarray
        Array of floats corresponding to the potentials representing the EEDF 
        in normalized probabilities.

    Notes
    -----
    The Druyvesteyn method requires the second derivative of the probe I-V 
    characteristic, which inherently amplifies noise and measurement errors.
    Therefore it is advisable to smooth the I-V prior to the use of this
    function.
    """
    
    probe_characteristic.check_validity()
    
    probe_characteristic.sort()
    
    characteristic_unique = probe_characteristic.get_unique()
    
    V_F = get_floating_potential(characteristic_unique)
    
    V_P = get_plasma_potential(characteristic_unique)
    
    # Obtain the correct EEDF energy range from the probe characteristic.
    _filter = (characteristic_unique.bias > V_F) & \
              (characteristic_unique.bias < V_P)
    energy = const.e * (V_P - characteristic_unique.bias[_filter])
    energy = energy.to(u.eV)
    
    # Below follows a very crude method of obtaining the second derivative of
    # the U-V curve. Any suggestions?
    dIdV = np.diff(characteristic_unique.current[_filter]) / \
           np.diff(characteristic_unique.bias[_filter])
    dIdV = np.append(dIdV, dIdV[-1])
    dIdV2 = np.diff(dIdV) / np.diff(characteristic_unique.bias[_filter])
    dIdV2 = np.append(dIdV2, dIdV2[-1])
        
    # Division by the Druyvesteyn factor. Since the result will be normalized
    # all constant values are irrelevant.
    probability = dIdV2 / energy
    
    # Integration of the EEDF for the purpose of normalization.
    integral = 0
    for i in np.arange(len(probability)):
        if(i > 0):
            integral = integral + probability[i] * np.abs(energy[i-1] - \
                                                          energy[i])    
    probability = probability / integral
    
    if(visualize):
        with quantity_support(): 
            plt.figure()
            plt.semilogy(energy, probability, c='k')
            plt.title("EEDF")
    
    return energy, probability