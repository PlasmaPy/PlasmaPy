__all__ = ["cold_plasma_function_solution"]

import numpy as np
import astropy.units as u
import astropy.constants as const
from plasmapy.formulary import parameters as pfp

from sympy import symbols, solve

def cold_plasma_function_solution(
    B: u.T,
    k: u.rad/u.m,
    omega_p: u.rad/u.s,
    omega_e:u.rad/u.s,
    theta: u.rad,
):
 
    #---#
    
    for arg_name in ("B"):
        value = locals()[arg_name].squeeze()
        if not (value.ndim == 0):
            raise ValueError(
                f"Argument '{arg_name}' must be a float or an integer."
                f"shape {value.shape}."
                )
        locals()[arg_name] = value
    
    k = k.squeeze()
    if not (k.ndim == 0  or k.ndim == 1):
        raise ValueError(
            f"Arguement 'k' needs to be a single value or a 1D array astropy Quantity,"
            f"got a value of shape {k.shpae}."
        )

    omega_p = omega_p.squeeze()
    omega_e = omega_e.squeeze()
    if not (omega_e.ndim == 0 or omega_e.ndim == 1):
        raise ValueError(
            f"Arguement 'omega_e' needs to be a single value or a single valued 1D array astropy Quantity,"   
            f"got value of shape {omega_e.shape}."
        )
        
    theta = theta.squeeze()
    theta = theta.to(u.radian)
    if not (theta.ndim == 0):
        raise ValueError(
            f"Arguement 'theta' needs to be a single value astropy Quantity,"
            f"got value of shape {theta.shape}."
        )
        
    #---#
    
    #---#
    
    k_dim = k.ndim
    if k_dim == 0:
        ck = np.tile(0*u.rad/u.s,1)
        ck[0] = k*const.c 
        k_int = True
    elif k_dim == 1:
        ck = np.tile(0*u.rad/u.s,len(k))
        for i in range(len(k)):
            ck[i] = k[i]*const.c 
        k_int = False
    else:
       raise ValueError(
            f"Arguement 'k' needs to be a single value or 1D array astropy Quantity,"
            f"got value of shape {k.shape}."
        )
       
    component_frequency = np.tile(0*u.rad/u.s, 2)
    component_frequency[0] = pfp.gyrofrequency(B=B, particle='H+', signed=False)
    component_frequency[1] = pfp.gyrofrequency(B=B, particle='e-', signed=True)
        
    lengths = []
    if omega_p.ndim == 0 and omega_e.ndim == 0:
        omega_int = True
        lengths.append(1)
        lengths.append(1)
    elif omega_p.ndim == 1 and omega_e.ndim == 1:
        omega_int = False
        lengths.append(len(omega_p))
        lengths.append(len(omega_e))
    else:
        raise ValueError(
            f"Arguement 'omega_p' and 'omega_e' need to be the same quantity type,"
            f"got value of shape {omega_p.shape} and {omega_e.shape}."
        )
    
    sum_len = min(lengths)
    
    plasma_proton = np.tile(0*u.rad/u.s, sum_len)
    plasma_electron = np.tile(0*u.rad/u.s, sum_len)
    
    if omega_int == False:
        for i in range(sum_len):
            plasma_proton[i] = omega_p[i]
            plasma_electron[i] = omega_e[i]
    elif omega_int == True:
        plasma_proton[0] = omega_p
        plasma_electron[0] = omega_e
    else:
        raise ValueError(
            f"Arguement 'omega_p' and 'omega_e' quantity type could not be deterimined,"
            f"got value of shape {omega_p.shape} and {omega_e.shape}."
        )
    
    #---#
    
    #---#
    
    w = 0*u.rad/u.s
    
    prot_w = (component_frequency[0])
    elec_w = (component_frequency[1])
    
    S = 1 - (plasma_proton[0]**2)/(w**2+(prot_w**2)) - (plasma_electron[0]**2)/(w**2+(elec_w**2))
    P = 1 - (plasma_proton[0]**2)/(w**2) - (plasma_electron[0]**2)/(w**2)
    D = ((prot_w)/(w))*((plasma_proton[0]**2)/(w**2-(prot_w**2))) + ((elec_w)/(w))*((plasma_electron[0]**2)/(w**2-(elec_w**2))) 
    
    print(S)
    
    R = S + D
    L = S - D
    
    A = S*(np.sin(theta)**2) + P*(np.cos(theta)**2)
    B = R*L*(np.sin(theta)**2) + P*S*(1 + np.cos(theta)**2)
    C = P*R*L
    
    print(A,B,C)

    #eq = A*((ck[0]/w)**4) - B*((ck[0]/w)**2)
    #print(eq)
    #sol = solve(eq, w)



    
    #print(sol)

    
    return 


inputs = {
   "B": 8.3e-9 * u.T,
   "k": 0.01 * u.rad / u.m,
   "omega_p": 1.6e6 * u.rad / u.s,
   "omega_e": 4.0e5 * u.rad / u.s,
   "theta": 30 * u.deg,
}

cold_plasma_function_solution(**inputs)

