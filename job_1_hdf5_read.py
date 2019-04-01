import numpy as np
import h5py 
f = h5py.File('mytestfile3.hdf5', 'r')
print(f.keys())
m=f['bigBOOL'][...]
f.close()