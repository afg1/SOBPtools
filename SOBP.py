"""
    use the equations from Jette & Chen, 2011 (Phys. Med. Biol. 56 N131) to produce the necessary ranges and weights to produce a SOBP using the pencil beam
    algorithm from Bortfeld. This includes their re-fudging of the value of p to generate a flat SOBP when simulated in Monte Carlo. 
    
    Version 0.0.1a: Working calculation of the SOBP, but I have had to use my own fudge factor to generate a flat SOBP as this code produces it (mainly to make a pretty picture). This only works for the specific calculation used below, and it is probably a better idea to use the lookup table in the function lookupP.
"""

import BortfeldBraggPeak as bfbp
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si

def lookupP(R0, chi, p0=1.77):
    alpha = 0.0022
    E = (R0/alpha)**(1/p0)
    print(E)

    Evals = np.array([50, 100, 150, 200, 250])
    p_15pc = np.array([1.48, 1.46, 1.43, 1.40, 1.34])
    p_20pc = np.array([1.45, 1.43, 1.40, 1.37, 1.32])
    p_25pc = np.array([1.43, 1.42, 1.39, 1.34, 1.29])
    p_30pc = np.array([1.43, 1.41, 1.37, 1.33, 1.27])
    p_35pc = np.array([1.42, 1.40, 1.36, 1.32, 1.26])
    p_40pc = np.array([1.41, 1.38, 1.35, 1.30, 1.24])
    if chi <= 0.15:
        interpolator = si.interp1d(Evals, p_15pc)
        return interpolator(E)
    elif chi <= 0.20:
        interpolator = si.interp1d(Evals, p_20pc)
        return interpolator(E)
    elif chi <= 0.25:
        interpolator = si.interp1d(Evals, p_25pc)
        return interpolator(E)
    elif chi <= 0.30:
        interpolator = si.interp1d(Evals, p_30pc)
        return interpolator(E)
    elif chi <= 0.35:
        interpolator = si.interp1d(Evals, p_35pc)
        return interpolator(E)
    elif chi <= 0.40:
        interpolator = si.interp1d(Evals, p_40pc)
        return interpolator(E)
    else:
        return 1.77 # the standard value, we can't really deal with energies this high yet






def calculateRanges(chi, R0, n):
    rk = []
    for k in np.arange(0, n+1, dtype=float):
        rk.append((1 - (( 1 - (k/n))*chi))*R0)
    return np.array(rk)


def calculateWeights(n, p):
    wk = []
    invp = 1.0/p
    power = 1.0 - invp
    for k in np.arange(0, n+1, dtype=float):
        if k == 0:
            wk.append(1.0 - ((1.0 - (1.0/(2.0*n)))**power ))
        elif k < n:
            wk.append(((1.0 - (1.0/n)*(k - 0.5))**power) - ((1 - (1.0/n)*(k + 0.5))**power))
        elif k == n:
            wk.append((1.0/(2.0*n))**power)
    return np.array(wk)


R0 = 15
chi = 0.15
lookupP(R0, chi)
p = 1.7 #lookupP(R0, chi) # Fudged this to get a nice flat SOBP with this model. It is unlikely to look like this when run through a Monte Carlo!
print(p)



rk = calculateRanges(chi, R0, 10)
wk = calculateWeights(10, p)
print(np.sum(wk))

print(rk)
print(wk)

wk *= 1000000000 *35

Zdata = np.linspace(0, 17, 1000)

DData = bfbp.GenerateBraggPeak(rk[0], wk[0], 0.20, 0.27, zlims=(0, 17))
plt.plot(Zdata, DData)
for r, w in zip(rk[1:], wk[1:]):
    DDatapart = bfbp.GenerateBraggPeak(r, w, 0.20, 0.27, zlims=(0,17))
    plt.plot(Zdata, DDatapart)
    DData += DDatapart


plt.plot(Zdata, DData)
plt.xlabel("Depth (cm)")
plt.ylabel("Dose (Gy)")
plt.savefig("SOBP_figure.pdf")
plt.show()