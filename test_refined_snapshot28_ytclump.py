"""

http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html

After getting the finely gridded cube from resample.py, use yt clump finder.


Last mod: 9 July 2018

NOTE
----
The yt clump finder was initially described in http://adsabs.harvard.edu/abs/2009ApJ...691..441S , but it's changed since then. What it does now is decompose into non-overlapping tiles (stored in a kd-tree), identify contours within a tile, and then connect them across tiles. It does this with an upper and lower bound on a field value, and looking for topologically connected sets.

With yt, connected sets can be identified. This enables the location and analysis of a hierarchy of clumps; in star formation, for instance, this allows the construction of diagrams describing the density at which fragmentation occurs.

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


def ytclumpfind(ds, dd, field=('density'), n_cut=100, step=10, N_cell_min=20, save=False, plot=True, saveplot=None):
    '''

    The way it's implemented now only works for single "density" field.

    Parameter
    ---------
    ds: yt StreamDataset

    dd: YTRegion

    field: tuple
        tuple of str
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
    if n_cut < 1:
        n_cut = 1.0    # to make sure whatever comes out after multiplicative by step won't be too small
    c_min = n_cut
    c_max = dd[field].max()

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
                prj.save('clumps0_' + str(int(step)) + '-' + vv + 'axis.png')
            else:
                prj.show()


    if plot:
        plotclumps(ds, saveplot=saveplot)


    return master_clump, leaf_clumps

if __name__ == "__main__"
  master10, leaf10 = ytclumpfind(ds, dd, ("density"),
                               n_cut=100,
                               step=10,
                               N_cell_min=20,
                               plot=True,
                               saveplot=True)

  print(master10.children)
  print(master10.children[0]['density'])    # children
  print(master10.children[0].children[0]['density'])   # sub-children

  # traverse the entire clump tree .
  for clump in master10:
    # print(clump.total_clumps)
    print(clump.clump_id)

  write_clump_index(master10, 0, "master10_clump_hierarchy.txt")
  write_clumps(master10, 0,  "master10_clumps.txt")

  import os
  os.system('cat *_clump_hierarchy.txt')
  os.system('cat *_clumps.txt')


  # look at leave prop.
  print("see more.. in http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html")

  for ind in range(len(leaf10)):
    print(leaf10[ind]["density"])
    print(leaf10[ind].quantities.total_mass())
    print(leaf10[ind].quantities.center_of_mass())


  master20, leaf20 = ytclumpfind(ds, dd, ("density"),
                               n_cut=100,
                               step=20,
                               N_cell_min=20, plot=True, saveplot=True)
  write_clump_index(master20, 0, "master20_clump_hierarchy.txt")
  write_clumps(master20, 0,  "master20_clumps.txt")

  master30, leaf30 = ytclumpfind(ds, dd, ("density"),
                               n_cut=100,
                               step=30,
                               N_cell_min=20, plot=True, saveplot=True)
  write_clump_index(master30, 0, "master30_clump_hierarchy.txt")
  write_clumps(master30, 0,  "master30_clumps.txt")

  master70, leaf70 = ytclumpfind(ds, dd, ("density"),
                                n_cut=100,
                                step=70,
                                N_cell_min=20, plot=True, saveplot=True)
  write_clump_index(master70, 0, "master70_clump_hierarchy.txt")
  write_clumps(master70, 0,  "master70_clumps.txt")

  master100, leaf100 = ytclumpfind(ds, dd, ("density"),
                                 n_cut=100,
                                 step=100,
                                 N_cell_min=20, plot=True, saveplot=True)
  write_clump_index(master100, 0, "master100_clump_hierarchy.txt")
  write_clumps(master100, 0,  "master100_clumps.txt")


  master200, leaf200 = ytclumpfind(ds, dd, ("density"),
                                 n_cut=100,
                                 step=200,
                                 N_cell_min=20, plot=True, saveplot=True)
  write_clump_index(master200, 0, "master200_clump_hierarchy.txt")
  write_clumps(master200, 0,  "master200_clumps.txt")


