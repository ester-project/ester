# We calculate the flatness with the Roche model
# calculate omk knowing omc and vice-versa

from numpy import *
from scipy.optimize import root


# we have to solve a cubic equation a-J2*a**3=1+J2+0.5*omk**2

def eps(omk):
	return omk**2/(2+omk**2)

def om_k(omc):
	khi=arcsin(omc)
	return sqrt(6*sin(khi/3)/omc-2)

omc=0.88
print 'omc=',omc,' omk=',om_k(omc)
