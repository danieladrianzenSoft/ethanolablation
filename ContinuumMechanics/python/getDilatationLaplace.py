#import math as m
import numpy as np
import mpmath as mp

def getDilatationLaplace(s,r,p):
	Rtot = p['Rtot']
	r0 = p['r0']
	K = p['K']
	lam = p['lambda']
	mu = p['mu']
	phi = p['phi']
	Q = p['Q']
	
	c = K * (lam + 2 * mu)
	D = (Q * phi) / (4 * mp.pi * (r0**2) * c * s)
	C = -(phi * r0 * s) / (4 * mu * K)
	
	al = ((s/c)**(1/2))
	#r = np.array(r,dtype=float)
	#r = mp.matrix(r)
	# vector = np.vectorize(np.float_)

	# e = (((D * (r0**2)) / r) * 
	# 	(np.sinh(al*(Rtot - r)) / ((1-(r0*C)) 
	# 	* (np.sinh(al*(Rtot-r0))
	# 	+ al*r0*np.cosh(al*(Rtot-r0))))))

	e = []
	for r_elem in r:
		new_e = (((D * (r0**2)) / r_elem) * 
			(mp.sinh(al*(Rtot - r_elem)) / ((1-(r0*C)) 
			* (mp.sinh(al*(Rtot-r0))
			+ al*r0*mp.cosh(al*(Rtot-r0))))))
		e.append(new_e)

	return e
