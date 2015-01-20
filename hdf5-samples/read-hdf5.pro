file = 'star.hdf5'
fid = h5f_open(file)
did = h5d_open(fid, '/star/T')
T = h5d_read(did)

print, 'Temperature at the center: ', T[0, 0]

h5d_close, did
h5f_close, fid
