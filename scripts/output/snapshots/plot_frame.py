import matplotlib.pyplot as plt
import h5py
snap =  h5py.File("000000000000001.h5",'r')

# show the state variables
snap.keys()

# example vorticity
plt.contourf(snap['q'][:])
plt.show()
