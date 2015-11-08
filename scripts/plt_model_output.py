import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt

simulation = "output"
pathi = simulation+"/snapshots/"
model =  h5py.File(simulation+"/setup.h5",'r')
diags =  h5py.File(simulation+"/diagnostics.h5",'r')


