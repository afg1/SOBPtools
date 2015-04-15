"""
My attempt to implement the Bragg Peak equation found in Bortfeld '97, Med Phys 24 (12) 2024-

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

#R0 = 13.5
#phi0 = 10000
#epsilon = 0.20
#sig = 0.27


def GenerateBraggPeak(R, phi, eps, sig, zlims=(0, 15), samples=1000):
    Dv = np.vectorize(D)
    zdata = np.linspace(zlims[0], zlims[1], samples)
    Ddata = Dv(R, phi, eps, sig, zdata)
    return Ddata
#Ddata2 = D_prec(R0, phi0, epsilon, sig, zdata)
#Ddata = Ddata/max(Ddata)


#Ddata, grad = sp.pbdv(-0.565, zdata)
#print(sp.pbdv(-0.565, 0.0))

#R0 = 12.5
#phi0 = 5000
#Ddata2 += D_prec(R0, phi0, epsilon, sig, zdata)

#plt.plot(zdata, Ddata)
##plt.plot(zdata, Ddata)
##plt.savefig("PristineBP.pdf")
#plt.show()