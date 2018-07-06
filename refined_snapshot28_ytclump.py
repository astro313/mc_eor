"""

http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html


After getting the finely gridded cube from resample.py, try yt clump finder again.

Last mod: 5 July 2018

NOTE
----
The yt clump finder was initially described in http://adsabs.harvard.edu/abs/2009ApJ...691..441S , but it's changed since then. What it does now is decompose into non-overlapping tiles (stored in a kd-tree), identify contours within a tile, and then connect them across tiles. It does this with an upper and lower bound on a field value, and looking for topologically connected sets.

With yt, connected sets can be identified. This enables the location and analysis of a hierarchy of clumps; in star formation, for instance, this allows the construction of diagrams describing the density at which fragmentation occurs.

"""

import numpy as np
import matplotlib.pyplot as plt

import yt
from yt.analysis_modules.level_sets.api import Clump, find_clumps, get_lowest_clumps
import h5py

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
field = ('density')   # can be a tuple of multiple fields
n_cut = 100

# multiplicative interval between contours
step = 50.0
N_cell_min = 20

# c_min = 10**np.floor(np.log10(dd[field]).min()  )
c_min = n_cut
c_max = 10**np.floor(np.log10(dd[field]).max()+1)
c_max = dd[field].max()

print "min/max value for finding contours: ", c_min, c_max

master_clump = Clump(dd, field)                       # this "base clump" just covers the whole domain.
master_clump.add_validator("min_cells", N_cell_min)   # weed out clumps with less than 20 cells.

# start finding
find_clumps(master_clump, c_min, c_max, step)


# Save the clump tree as a reloadable dataset
fn = master_clump.save_as_dataset(fields=["density"]),  # "particle_mass"])
# # To reload the clump dataset
# cds = yt.load(fn)
# leaf_clumps_reloaded = cds.leaves
import os
os.system('ls -lrt')

leaf_clumps = get_lowest_clumps(master_clump)          # traverse clump hierarchy to get list of all 'leaf' clumps, which are the individual clumps that have no children of their own
print(leaf_clumps[0]["density"])
print(leaf_clumps[0].quantities.total_mass())
print(leaf_clumps[0].quantities.center_of_mass())


# the master clump will represent the top of a hierarchy of clumps
print(master_clump.children[0]['density'])
for clump in master_clump:
    print(clump.clump_id)

print("see more.. in http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html")


# write function to overplot the clumps found (specifically leaf_clumps)
prj = yt.ProjectionPlot(ds, 0, ("density"),
                        center='c')
prj.annotate_clumps(leaf_clumps)
prj.save('clumps0-xaxis.png')
# prj.show()



# ----------------------------------------- will look at after meeting .......
# multiplicative interval between contours
step = 50.0
N_cell_min = 20

# c_min = 10**np.floor(np.log10(dd[field]).min()  )
c_min = n_cut
c_max = 10**np.floor(np.log10(dd[field]).max()+1)
c_max = dd[field].max()

print "min/max value for finding contours: ", c_min, c_max

master_clump = Clump(dd, field)                       # this "base clump" just covers the whole domain.
master_clump.add_validator("min_cells", N_cell_min)   # weed out clumps with less than 20 cells.

# start finding
find_clumps(master_clump, c_min, c_max, step)


# Save the clump tree as a reloadable dataset
fn = master_clump.save_as_dataset(fields=["density"]),  # "particle_mass"])
# # To reload the clump dataset
# cds = yt.load(fn)
# leaf_clumps_reloaded = cds.leaves
import os
os.system('ls -lrt')

leaf_clumps = get_lowest_clumps(master_clump)          # traverse clump hierarchy to get list of all 'leaf' clumps, which are the individual clumps that have no children of their own
print(leaf_clumps[0]["density"])
print(leaf_clumps[0].quantities.total_mass())
print(leaf_clumps[0].quantities.center_of_mass())


# the master clump will represent the top of a hierarchy of clumps
print(master_clump.children[0]['density'])
for clump in master_clump:
    print(clump.clump_id)

print("see more.. in http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html")


# write function to overplot the clumps found (specifically leaf_clumps)
prj = yt.ProjectionPlot(ds, 0, ("density"),
                        center='c')
prj.annotate_clumps(leaf_clumps)
prj.save('clumps0-xaxis.png')
# prj.show()


# --- Save cloud properties ---
cloud_keys = "...."
