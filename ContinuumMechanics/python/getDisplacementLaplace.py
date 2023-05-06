import math as m

def getDisplacementLaplace(s, r, p):
	Rtot = p['Rtot']
	r0 = p['r0']
	K = p['K']
	lam = p['lambda']
	mu = p['mu']
	phi = p['phi']
	Q = p['Q']
	
	c = K * (lam + 2 * mu)
	al = (s/c)**(1/2)
	
	D = (Q*phi)/(4*m.pi*(r0**2)*c*s)
	C = -(phi*r0*s)/(4*mu*K)
	F = ((D*r0^2)/(((1-(r0*C))*(m.sinh(al*(Rtot-r0))))
		+al*r0*m.cosh(al*(Rtot-r0))))
	A = ((r0/al)*m.cosh(al*(r0-Rtot))-((1/(al**2))+((2*mu+lam)/(4*mu))*(r0**2))*m.sinh(al*(r0-Rtot)))
	
	u = F*(-m.sinh(al*(Rtot-r))/(al^2*r^2))-F*m.cosh(al*(Rtot-r))/(al*r)+(F*A)/(r**2)

	return u
