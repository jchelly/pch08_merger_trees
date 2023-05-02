#!/bin/env python

import numpy as np
import h5py
import matplotlib.pyplot as plt

def read_progenitors(filename):
    
    prog_t = np.dtype([
        ("tree_index", np.int32),
        ("mphalo", np.float32),
        ("mprog", np.float32),
    ])
    data = np.loadtxt(filename, dtype=prog_t)
    return data

def make_plot():

    eps = read_progenitors("../output/eps_mphalo1e12_z4.txt")
    pch = read_progenitors("../output/pch_mphalo1e12_z4.txt")
    ajb = read_progenitors("../output/ajb2022_mphalo1e12_z4.txt")
    galform = read_progenitors("../output/galform_benson2022.txt")

    for data, ls, label in ((eps, "r-", "EPS"),
                            (pch, "g-", "PCH"),
                            (ajb, "b-", "Benson 2022"),
                            (galform, "k-", "Galform (Benson 2022)")):
        tree_index = data["tree_index"]
        mphalo = data["mphalo"]
        mprog = data["mprog"]
        f = mprog/mphalo
        bins = np.logspace(-5.0, 0.0, 50)
        plt.hist(f, bins=bins, label=label, histtype="step")

    plt.xlabel(r"$M_{prog}/M_{halo}$")
    plt.xscale("log")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.legend(loc="upper right")
    plt.title(r"100 halos of $10^{12} h^{-1} M_{\odot}$, halo at z=0, progenitors at z=4")

if __name__ == "__main__":
    make_plot()
    plt.show()
