import numpy as np
from getDilatationLaplace import getDilatationLaplace
from getDisplacementLaplace import getDisplacementLaplace
import matplotlib.pyplot as plt
import mpmath as mp

K = 3.8e-11
r0 = 0.035
Rtot = 1
Q = (0.1/1e3)/60
lam = 13.16*1e4
mu = 6.58*1e4
phi = 0.2
numr = 20
R = np.linspace(r0,Rtot,numr)
r_mass = np.concatenate((np.linspace(0,r0-1e-3,20),(np.linspace(r0,Rtot,numr))))
T = np.linspace(0.01,1,6*60)

p = {'K': K, 'r0': r0, 'Rtot': Rtot,
	'Q': Q, 'lambda': lam, 'mu': mu,
	'phi': phi}

#el = lambda s: getDilatationLaplace(s,r,p)
ep = lambda t: mp.invertlaplace(el,t,
                                      method='dehoog',
                                      dps=5, 
                                      degree=3)

e = []
for r in R:
	for t in T:
		el = lambda s: getDilatationLaplace(s,r,p)
		e_temp = mp.invertlaplace(el, t, method='dehoog',
                                	dps=5, 
                                    degree=3)
		e.append(e_temp)
#e = np.array([[ep(t,r) for t in T] for r in R])
#up = list()

#for i in range(0,len(r)):


# b = [2,6]
# a = np.array([1,4,5,7,3])
# res = np.searchsorted(a,6)

print(e)