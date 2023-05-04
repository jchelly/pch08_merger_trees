#!/bin/env python

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os.path


def read_progenitors(filename):

    txt_filename = filename
    hdf5_filename = filename+".hdf5"

    # Check if HDF5 version of file is newer
    txt_mtime = os.path.getmtime(txt_filename)
    try:
        hdf5_mtime = os.path.getmtime(hdf5_filename)
    except FileNotFoundError:
        have_hdf5 = False
    else:
        have_hdf5 = hdf5_mtime > txt_mtime

    prog_t = np.dtype([
        ("tree_index", np.int32),
        ("mphalo", np.float32),
        ("mprog", np.float32),
        ("jlevel", np.int32),
        ("z", np.float32),
    ])

    if have_hdf5:
        print(f"Using {hdf5_filename}")
        # Get cached HDF5 version of data
        with h5py.File(hdf5_filename, "r") as infile:
            data = infile["progenitors"][...]
    else:
        print(f"Need to regenerate HDF5 file {hdf5_filename}")
        # Read the input file and make a new HDF5 file
        data = np.loadtxt(txt_filename, dtype=prog_t)
        with h5py.File(hdf5_filename, "w") as outfile:
            outfile["progenitors"] = data

    # Sort by redshift
    data = np.sort(data, order="z")

    # Find unique redshifts
    unique_z, unique_offset, unique_count = np.unique(data["z"], return_index=True, return_counts=True)
    
    # Construct dictionary for each redshift
    progenitors = {}
    for z, offset, count in zip(unique_z, unique_offset, unique_count):
        progenitors[z] = data[offset:offset+count]
        assert np.all(progenitors[z]["z"] == z)

    return progenitors


def round_z(prog):
    prog_new = {}
    for z in prog:
        z_rounded = np.around(z)
        prog_new[z_rounded] = prog[z]
    return prog_new


def make_plot():

    # Read progenitor info
    trees_ajb_prog = read_progenitors("/cosma7/data/dp004/jch/benson2022_m1e12_10000.dat")
    trees_pch_prog = read_progenitors("/cosma7/data/dp004/jch/pch2008_m1e12_10000.dat")
    galform_ajb_prog = read_progenitors("/cosma7/data/dp004/jch/galform_benson2022_m1e12_10000.dat")

    # Round off redshifts
    trees_ajb_prog = round_z(trees_ajb_prog)
    trees_pch_prog = round_z(trees_pch_prog)
    galform_ajb_prog = round_z(galform_ajb_prog)

    prog_z = sorted([p for p in list(trees_ajb_prog) if p != 0.0])

    # Make a plot for each redshift
    for i, z in enumerate(prog_z):

        plt.subplot(2,2,i+1)

        for label, prog in (("Benson 2022 (modified code)", trees_ajb_prog),
                            ("PCH08 (original code)", trees_pch_prog),
                            ("Benson 2022 (Galform)", galform_ajb_prog)):
            mphalo = prog[z]["mphalo"]
            mprog = prog[z]["mprog"]
            f = mprog/mphalo
            bins = np.logspace(-5.0, 0.0, 50)

            num_per_bin, edges = np.histogram(f, bins=bins)
            centres = np.sqrt(edges[1:]*edges[:-1])
            errors = np.sqrt(num_per_bin)
            plt.errorbar(centres, num_per_bin, yerr=errors, label=label)
            #plt.hist(f, bins=bins, histtype="step", label=label)

        plt.xlabel(r"$\rm{}M_{prog}/M_{halo}$")
        plt.xscale("log")
        plt.ylabel("Frequency")
        plt.yscale("log")
        plt.ylim(1.0e2, 1.0e6)
        plt.title(f"z = {z}")
        if i == 0:
            plt.legend(loc="lower left")

if __name__ == "__main__":

    plt.figure(figsize=(8,8))
    make_plot()
    plt.tight_layout()
    plt.savefig("progenitor_masses.pdf")
    plt.show()
