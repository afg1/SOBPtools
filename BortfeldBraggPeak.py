"""
My attempt to implement the Bragg Peak equation found in Bortfeld '97, Med Phys 24 (12) 2024-

Version 0.0.1a: Working implementation, by default we use the 'inaccurate' version of the dose calculation, mainly for speed. May add another call to enable use of the accurate version if required. Note: these equations are valid for water only, and return the dose in Gray via the conversion factor defined globally (calculations are done in MeV and cm).

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

toGray = 1.602E-10

def D(R0, phi0, epsilon, sig, z):
    """
    This is the very specialised equation for water only. It's equation 28/29 in the paper
    """
    if z < (R0 - 10*sig):
        fac = (phi0)/(1.0 + 0.012*R0)
        term1 = 17.93*((R0 - z)**-0.435)
        term2 = ((0.444 + 31.7*epsilon)/R0)*((R0-z)**0.565)
        return fac*(term1 + term2) * toGray
    
    elif  z < (R0 + 5*sig):
        D565, grad565 = sp.pbdv(-0.565, -((R0 - z)/sig))
        D1565, grad1565 = sp.pbdv(-1.565, -((R0 - z)/sig))
        
        frontfac = ((np.exp( (-(R0 - z)**2)/(4.0*(sig**2))) * (sig**0.565))/(1.0 + 0.012*R0))*phi0
        bracfac = 11.26*D565/sig + ((0.157 + 11.26*epsilon)/R0)*D1565
        return frontfac * bracfac * toGray
    
    else:
        return 0.0

def D_prec(R0, phi0, epsilon, sig, z):
    D565, grad565 = sp.pbdv(-0.565, -((R0 - z)/sig))
    D1565, grad1565 = sp.pbdv(-1.565, -((R0 - z)/sig))
    
    frontfac = ((np.exp( (-(R0 - z)**2)/(4.0*(sig**2))) * (sig**0.565))/(1.0 + 0.012*R0))*phi0
    bracfac = 11.26*D565/sig + ((0.157 + 11.26*epsilon)/R0)*D1565
    return frontfac * bracfac * toGray




def GenerateBraggPeak(R, phi, eps, sig, zlims=(0, 15), samples=1000):
    Dv = np.vectorize(D)
    zdata = np.linspace(zlims[0], zlims[1], samples)
    Ddata = Dv(R, phi, eps, sig, zdata)
    return Ddata
