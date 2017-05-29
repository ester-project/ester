# radial integrator : compute the primitive that vanishes at the origin

import numpy as np

def integre(f,a):

# computes the jfirst grid points
    jfirst=np.zeros(a.ndomains,dtype=np.int)

    for i in np.arange(a.ndomains-1)+1:
       jfirst[i]=jfirst[i-1]+a.npts[i]

# fill the matrix with the derivation-matrix
    A=a.D
# Impose the continuity of the primitive at each domain junction
    for idom in np.arange(a.ndomains-1)+1:
        if idom >= 1:
           A[jfirst[idom],:]=0.0
           A[jfirst[idom],jfirst[idom]-1]=-1.0
           A[jfirst[idom],jfirst[idom]]=1.0
           f[jfirst[idom]]=0.0

# Impose that the primitive vanishes at the origin
    A[0,:]=0.
    A[0,0]=1.
    f[0]=0.

# Solve for the RHS
    x = np.linalg.solve(A, f)

    return x

# and get the primitive on all the Gauss-Lobatto grid(s)
