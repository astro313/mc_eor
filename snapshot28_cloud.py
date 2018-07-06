"""

Last mod: 5 July 2018

History:
5 July 2018:
    - modified Juan's snippet to test snapshot 28, using "Density" field

see also https://github.com/astro313/mc_eor/blob/master/test_yt_clump.ipynb

"""

from yt.config import ytcfg
ytcfg["yt", "__withinreason"] = "True"
import os
#import Pyro4
import uuid

from yt.mods import *
import numpy as np
import math
import pylab as P
import matplotlib.pyplot as plt

import yt
import h5py

kpc = 3.0856e21  # cm
pc = 3.0856e18  # cm
km = 1.0e5      # cm
Myr = 3.1556e13  # s
mp = 1.6726e-24  # g
mu = 1.2924      # what is this?
kb = 1.3806e-16  # erg K-1
GNewton = 6.6743e-8   # cm3 g-1 s-2
Msun = 1.9884e33   # g
mm = mu * mp

Plot_stuff = False
debug = False

# convert from code unit density to g/cc (depending on how fetch_gal.py is implemented.)
convert_unit = True

if debug:
    namefile = "output/output_00028/info_00028.txt"
    myfile = namefile

    print "Loading file,", myfile
    pf = load(myfile)


f = h5py.File("snapshot28_center_densityfield_resampled.h5", "r")
density = f["density"].value

if convert_unit:
    import pymses
    from pymses.utils import constants as C

    ro = pymses.RamsesOutput("output", 28)
    factor = ro.info["unit_density"].express(C.H_cc)
    density *= factor
    print density.max()

data = dict(density=density)
ds = yt.load_uniform_grid(data, f["density"].shape)

dd = ds.all_data()
# dd['density']

# Define density threshold to cut the clouds.
field = 'density'
n_cut = 100
rho_cut = n_cut

c_max = dd[field].max()

# Main call of this script.
# This function extract_connected_sets,
# generates a tree of connected structures w/in those structures.
# Here I only request 1 structure, which is the
# entire cloud.

cons, contours = dd.extract_connected_sets(field, 1, rho_cut, c_max)

# Keep the clouds that are at least N_cell_min cells in volume and drop the rest.
N_cell_min = 20
num_contours = len(contours[0])
fake_clouds = 0
all_clouds = np.ones(num_contours)
cloud_list = []

for i in range(num_contours):
    obj = contours[0][i]['CellMass']

# ---------------------------------------------------------------------------
# KeyError                                  Traceback (most recent call last)
# <ipython-input-58-db901a3a4ee2> in <module>()
# ----> 1 contours[0][0]['CellMass']

# KeyError: 0
    if obj.size < N_cell_min:
        fake_clouds += 1
        all_clouds[i] = -1
    else:
        obj = contours[0][i]
        cloud_list.append(obj)

# number of real clouds
num_real_clouds = num_contours - fake_clouds

print " Total number of contours found       : ", num_contours
print " The number of unresolved clouds is   : ", fake_clouds
print " The total number of clouds is        : ", num_real_clouds




# --- Save cloud properties ---
cloud_keys = "...."




# ---