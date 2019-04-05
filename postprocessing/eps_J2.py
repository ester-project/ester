# We calculate the flatness with the Roche model taking into account the quadrupolar
# moment of the gravitational potential

from numpy import *
from scipy.optimize import root


# we have to solve a cubic equation a-J2*a**3=1+J2+0.5*omk**2

def a(x,J2,omk):
	return x-J2*x**3-1.-J2-0.5*omk**2

J2=0.002108
omk=0.8366
print a(0.9,J2,omk)
rt= root(a, 1., args=(J2,omk )).x[0]
rt0= root(a, 1., args=(0.,omk )).x[0]

eps=1.-1./rt
eps0=1.-1./rt0
print 'eps = ',eps,'eps (J2=0) = ',eps0
print 'J2 effect',(eps-eps0)/eps0
