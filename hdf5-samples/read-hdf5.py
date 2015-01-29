import h5py

f = h5py.File('star.hdf5', 'r')  # open the file 'star.hdf5' ('r' = read only)
T = f['/star/T'][:]              # read the temperature field
n = f['/star'].attrs['ndomains'] # read the number of domains

print "T (center): " + str(T[0][0])
print "T (equator): " + str(T[0][-1])
