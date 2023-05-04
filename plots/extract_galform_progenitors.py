#!/bin/env python

import h5py
import numpy as np
import virgo.util.match as match

filename="/cosma7/data/dp004/jch/galform_output/benson2022/galaxies.hdf5"
with h5py.File(filename, "r") as infile:
    mhalo = infile["Halo_Trees/mhalo"][...]
    snap_num = infile["Halo_Trees/SnapNum"][...]
    out_snap = infile["Output_Times/snapnum_out"][...]
    out_z    = infile["Output_Times/zout"][...]

# Find redshift of each halo
ptr = match.match(snap_num, out_snap)
assert np.all(ptr>=0)
halo_z = out_z[ptr]

# Write output
for i in range(len(mhalo)):
    print(f"0 1.0000000E+12 {mhalo[i]:.6e} {snap_num[i]} {halo_z[i]}")
