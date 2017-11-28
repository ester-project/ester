import os
import argparse
import	sys
import math
import numpy as np
import os.path

parser = argparse.ArgumentParser()

parser.add_argument('filename')
args = parser.parse_args()
filename=args.filename

gen=np.loadtxt(filename)

N=len(gen[:,0])
R=np.ones(N)
m=np.ones(N)
Z=0.02*np.ones(N)
T=np.zeros(N)
rho=np.zeros(N)
X=np.zeros(N)
kappa_gen=np.zeros(N)
kappa_est=np.zeros(N)

if os.path.isfile('output_temp.txt'):

	os.system('rm output_temp.txt')


os.system ('touch output_temp.txt') # Fichier d'output de pme.cpp
os.system('make')


for k in range(0,N):

	print "iteration ", k + 1, " / ", N
	T[k]=np.power(10,gen[k,3])
	rho[k]=np.power(10,gen[k,4])
	X[k]=gen[k,19]
	kappa_gen[k]=np.power(10,gen[k,13])
	R[k]=np.power(10,gen[k,1])
	m[k]=gen[k,2]


	os.system('./pme ' +'%10.9e' % rho[k] + ' %10.9e' % X[k] + ' %10.9e' % T[k])

est=np.loadtxt('output_temp.txt')
os.system('rm output_temp.txt') 

for i in range(0,len(est)): # Recup kappa_est depuis output_temp.txt

	kappa_est[i]=est[i]

np.savetxt('output.txt',(np.column_stack((m/(5.*1.98855e33),R,X,Z,T,rho,kappa_gen, kappa_est))))



