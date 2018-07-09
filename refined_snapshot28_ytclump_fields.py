"""

http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html

After getting the finely gridded cube from resample_fields.py, use yt clump finder.


Last mod: 9 July 2018

NOTE
----
The yt clump finder was initially described in http://adsabs.harvard.edu/abs/2009ApJ...691..441S , but it's changed since then. What it does now is decompose into non-overlapping tiles (stored in a kd-tree), identify contours within a tile, and then connect them across tiles. It does this with an upper and lower bound on a field value, and looking for topologically connected sets.

With yt, connected sets can be identified. This enables the location and analysis of a hierarchy of clumps; in star formation, for instance, this allows the construction of diagrams describing the density at which fragmentation occurs.

- extract_connected_sets() differ slightly from Clump()
    - extract_connected_sets(): get topologically connected sets. These sets are identified by examining cells between two threshold values and connecting them.
    - Clump(): uses a contouring algorithm to identified topologically disconnected structures within a dataset. This works by first creating a single contour over the full range of the contouring field, then continually increasing the lower value of the contour until it reaches the maximum value of the field. As disconnected structures are identified as separate contours, the routine continues recursively through each object, creating a hierarchy of clumps. Individual clumps can be kept or removed from the hierarchy based on the result of user-specified functions, such as checking for gravitational boundedness.

"""

import numpy as np
import matplotlib.pyplot as plt

import yt
from yt.analysis_modules.level_sets.api import Clump, find_clumps, get_lowest_clumps, write_clump_index, write_clumps
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


f = h5py.File("snapshot28_center_fields012345-15_resampled.h5", "r")
density = f["rho"].value      # careful sometimes i used "density" (e.g., resample.py), see resample_fields.py to make sure
H2 = f["H2"].value
Pressure = f["P"].value
# vel = f["vel"].value   # to implement in resample_fields.py

if convert_unit:
    import pymses
    from pymses.utils import constants as C

    ro = pymses.RamsesOutput("output", 28)
    factor = ro.info["unit_density"].express(C.H_cc)
    density *= factor
    print density.max()

data = dict(density=density, H2=H2)
ds = yt.load_uniform_grid(data, f["rho"].shape)
dd = ds.all_data()


def ytclumpfind_H2(ds, dd, field, n_cut=100, step=10, N_cell_min=20, save=False, plot=True, saveplot=None):
    '''

    The way it's implemented now only works for single "density" field.

    Parameter
    ---------
    ds: yt StreamDataset

    dd: YTRegion

    field: tuple
        tuple of str, e.g., ("gas", "density") or ("io", "....") or ("gas", "averaged_density"), etc... or just ("density") or ("stream", "density")

    n_cut: float or int
        defines lowest contour level to start searching from.
        e.g., density to cut the clouds.

    step: int
        multiplicative interval between subsequent contours

    N_cell_min: int
        min. number of cell s.t. clump identified is not spurious

    save: Boolean
        if true, save the clump tree as a reloadable dataset

    plot: Boolean
        if true, will plot leaf clump

    saveplot: boolean
        if true, will save figure instead of showing it


    Return
    ------
    master_clump: the top of a hierarchy of clumps
    leaf_clumps: list of individual clumps that have no children of their own


    '''

    if plot:
        assert saveplot is not None

    # c_min = 10**np.floor(np.log10(dd[field]).min()  )
    # c_max = 10**np.floor(np.log10(dd[field]).max()+1)
    if n_cut < 1.e-5:
        n_cut = 1.0    # to make sure whatever comes out after multiplicative by step won't be too small
    c_min = n_cut
    c_max = (dd[field]).max()

    # assert np.isnan(c_min) is False and np.isnan(c_max) is False

    print "min/max value for finding contours: ", c_min, c_max

    master_clump = Clump(dd, field)                       # this "base clump" just  covers the whole domain.
    master_clump.add_validator("min_cells", N_cell_min)   # weed out clumps < N_cell_min cells.

    find_clumps(master_clump, c_min, c_max, step)

    if save:

        fn = master_clump.save_as_dataset(fields=list(field)),  # "particle_mass"])
        # # To reload the clump dataset
        # cds = yt.load(fn)
        # leaf_clumps_reloaded = cds.leaves

    # traverse clump hierarchy to get list of all 'leaf' clumps, which are the individual clumps that have no children of their own
    leaf_clumps = get_lowest_clumps(master_clump)


    def plotclumps(ds, field=field, saveplot=saveplot):
        """ overplot the clumps found (specifically the leaf_clumps) along 3 images, each created by projecting onto x-, y-, and z-axis. """

        axes = {'0': 'x', '1': 'y', '2': 'z'}

        for kk, vv in axes.iteritems():

            prj = yt.ProjectionPlot(ds,
                                    int(kk),
                                    field,
                                    center='c')
            prj.annotate_clumps(leaf_clumps)
            if saveplot:
                prj.save('clumps1_' + str(int(step)) + '-' + vv + 'axis.png')
            else:
                prj.show()


    if plot:
        plotclumps(ds, saveplot=saveplot)


    return master_clump, leaf_clumps


# -------------- decide on n_cut --------------
# may be useful to plot the 3D: http://yt-project.org/doc/visualizing/volume_rendering.html, look at transfer functions

import yt
import numpy as np
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.api import Scene, VolumeSource


# Or see fig 6 Pallottini 2017 for n_H2 cut as starting point, lower right panel

# in units of nH2/cc
n_cut_1 = 10**0.5
n_cut_2 = 10**-1.5



# -------------- run clump finder -------------


master10, leaf10 = ytclumpfind_H2(ds, dd, ("density"),
                               n_cut=n_cut_1,
                               step=10,
                               N_cell_min=20,
                               plot=True,
                               saveplot=True)

import pdb; pdb.set_trace()

print(master10.children)
print(master10.children[0]['density'] * master10.children[0]['H2'])    # children
print(master10.children[0].children[0]['density'] * master10.children[0].children[0]['H2'])   # sub-children

# traverse the entire clump tree .
for clump in master10:
    # print(clump.total_clumps)
    print(clump.clump_id)

write_clump_index(master10, 0, "master10_clump_hierarchy_H2.txt")
write_clumps(master10, 0,  "master10_clumps_H2.txt")

import os
os.system('cat *_clump_hierarchy_H2.txt')
os.system('cat *_clumps_H2.txt')


for ind in range(len(leaf10)):
    print(leaf10[ind]["density"] * leaf10[ind]["H2"])
    print(leaf10[ind].quantities.total_mass())
    print(leaf10[ind].quantities.center_of_mass())
    print(leaf10[ind]["clump"])
    print(leaf10[ind]["grid"])


master20, leaf20 = ytclumpfind_H2(ds, dd, ("density"),
                               n_cut=n_cut_1,
                               step=20,
                               N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master20, 0, "master20_clump_hierarchy_H2.txt")
write_clumps(master20, 0,  "master20_clumps_H2.txt")

master30, leaf30 = ytclumpfind_H2(ds, dd, ("density"),
                               n_cut=n_cut_1,
                               step=30,
                               N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master30, 0, "master30_clump_hierarchy_H2.txt")
write_clumps(master30, 0,  "master30_clumps_H2.txt")

master70, leaf70 = ytclumpfind_H2(ds, dd, ("density"),
                                n_cut=n_cut_1,
                                step=70,
                                N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master70, 0, "master70_clump_hierarchy_H2.txt")
write_clumps(master70, 0,  "master70_clumps_H2.txt")

master100, leaf100 = ytclumpfind_H2(ds, dd, ("density"),
                                 n_cut=n_cut_1,
                                 step=100,
                                 N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master100, 0, "master100_clump_hierarchy_H2.txt")
write_clumps(master100, 0,  "master100_clumps_H2.txt")


master200, leaf200 = ytclumpfind_H2(ds, dd, ("density"),
                                 n_cut=n_cut_1,
                                 step=200,
                                 N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master200, 0, "master200_clump_hierarchy_H2.txt")
write_clumps(master200, 0,  "master200_clumps_H2.txt")



# --- repeat for n_cut_2 ---

master10, leaf10 = ytclumpfind_H2(ds, dd, ("density"),
                               n_cut=n_cut_2,
                               step=10,
                               N_cell_min=20,
                               plot=True,
                               saveplot=True)
write_clump_index(master10, 0, "master10_clump_hierarchy_H2.txt")
write_clumps(master10, 0,  "master10_clumps_H2.txt")

master20, leaf20 = ytclumpfind_H2(ds, dd, ("density"),
                               n_cut=n_cut_2,
                               step=20,
                               N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master20, 0, "master20_clump_hierarchy_H2.txt")
write_clumps(master20, 0,  "master20_clumps_H2.txt")

master30, leaf30 = ytclumpfind_H2(ds, dd, ("density"),
                               n_cut=n_cut_2,
                               step=30,
                               N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master30, 0, "master30_clump_hierarchy_H2.txt")
write_clumps(master30, 0,  "master30_clumps_H2.txt")

master70, leaf70 = ytclumpfind_H2(ds, dd, ("density"),
                                n_cut=n_cut_2,
                                step=70,
                                N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master70, 0, "master70_clump_hierarchy_H2.txt")
write_clumps(master70, 0,  "master70_clumps_H2.txt")

master100, leaf100 = ytclumpfind_H2(ds, dd, ("density"),
                                 n_cut=n_cut_2,
                                 step=100,
                                 N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master100, 0, "master100_clump_hierarchy_H2.txt")
write_clumps(master100, 0,  "master100_clumps_H2.txt")


master200, leaf200 = ytclumpfind_H2(ds, dd, ("density"),
                                 n_cut=n_cut_2,
                                 step=200,
                                 N_cell_min=20, plot=True, saveplot=True)
write_clump_index(master200, 0, "master200_clump_hierarchy_H2.txt")
write_clumps(master200, 0,  "master200_clumps_H2.txt")
