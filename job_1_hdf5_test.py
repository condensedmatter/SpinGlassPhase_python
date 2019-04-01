import h5py
import numpy

'''
f = h5py.File("mytestfile3.hdf5", "a")
 
s=numpy.random.random(100)

dset = f.create_dataset("rand",data=s)



dset[:]=s


f.close()

'''


d1 = np.random.random(size = (1000,20))
d2 = np.random.random(size = (1000,200))
s=np.random.choice(a=[False, True], size=300, p=[0.5,0.5])

hf = h5py.File('data.h5', 'w')

hf.create_dataset('dataset_1', data=d1)
da2=hf.create_dataset('dataset_2', data=d2)
dataset_s=hf.create_dataset('dataset_s', data=s)

print(da2[:])

#hf.close()
