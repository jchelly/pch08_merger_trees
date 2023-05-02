#!/bin/env python

import h5py
import numpy as np

filename="/cosma/home/jch/Codes/Galform/git/galform/output/benson2022/galaxies.hdf5"
with h5py.File(filename, "r") as infile:
    mhalo = infile["Halo_Trees/mhalo"][...]
    snapnum = infile["Halo_Trees/SnapNum"][...]

for i in range(len(mhalo)):
    if snapnum[i] == 61:
        print("0  1.0e12  ", mhalo[i])
