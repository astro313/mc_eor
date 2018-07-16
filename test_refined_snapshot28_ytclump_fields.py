"""

Potential BUG/issue/inconsisntency (may need to update code):  density should be multipled by factor = get_units(ro=ro)['rho'][0]      # 1/cm^3 (not H/cm^3) if we use Fig 6 of Pallottini+17 as n_cut.


http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html

After getting the finely gridded cube from resample_fields.py, use yt clump finder.


Last mod: 16 July 2018

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

fold_out = 'test_png_28/'


# convert from code unit density to g/cc (depending on how fetch_gal.py is
# implemented.)
convert_unit = True

f = h5py.File("snapshot28_center_fields0123456-15_resampled.h5", "r")
# careful sometimes i used "density" (e.g., resample.py), see
# resample_fields.py to make sure
density = f["rho"].value
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

# make a derived field, call h2density (for yt Clump() to work)
def _h2density(field, data):
    try:
        return data["density"] * data["H2"]
    except:
        return data[("stream", "density")] * data[("stream", "H2")]

print dd['H2'].max()
print dd['density'].max()
print (dd['H2'] * dd['density']).max()

from yt.units import dimensions
# The global yt.add_field() function is for adding a field for every
# subsequent dataset that is loaded in a particular python session,
# whereas ds.add_field() (add_field()) will only add it to dataset ds.
ds.add_field(("stream", "h2density"), function=_h2density, units="g/cm**3")
print dd['h2density'].max()

assert (dd['H2'] * dd['density']).max() == dd['h2density'].max()


def ytclumpfind_H2(ds, dd, field, n_cut, c_max=None, step=3, N_cell_min=10, save=False, plot=True, saveplot=None, fold_out='./'):
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

    fold_out: str
        directory to save figures


    Return
    ------
    master_clump: the top of a hierarchy of clumps
    leaf_clumps: list of individual clumps that have no children of their own


    '''

    if plot:
        assert saveplot is not None

    # c_min = 10**np.floor(np.log10(dd[field]).min()  )
    # c_max = 10**np.floor(np.log10(dd[field]).max()+1)
    # if n_cut < 1.e-5:
    #     n_cut = 1.0    # to make sure whatever comes out after multiplicative by step won't be too small

    if step <= 1.0:
        step += 1
    c_min = n_cut

    if c_max is None:
        c_max = (dd[field]).max()

    # assert np.isnan(c_min) is False and np.isnan(c_max) is False
    print "min/max value for finding contours: ", c_min, c_max

    # this "base clump" just  covers the whole domain.
    master_clump = Clump(dd, field)
    # weed out clumps < N_cell_min cells.
    master_clump.add_validator("min_cells", N_cell_min)

    find_clumps(master_clump, c_min, c_max, step)

    if save:

        fn = master_clump.save_as_dataset(
            fields=list(field)),  # "particle_mass"])
        # # To reload the clump dataset
        # cds = yt.load(fn)
        # leaf_clumps_reloaded = cds.leaves

    # traverse clump hierarchy to get list of all 'leaf' clumps, which are the
    # individual clumps that have no children of their own
    leaf_clumps = get_lowest_clumps(master_clump)

    # print "printing sth here...", ds.tree["clump", field]
    # # *** AttributeError: 'StreamDataset' object has no attribute 'tree'

    def plotclumps(ds, leaf_clumps, master_clump, field=field, saveplot=saveplot, fold_out=fold_out):
        """ overplot the clumps found (specifically the leaf_clumps) along 3 images, each created by projecting onto x-, y-, and z-axis. """

        axes = {'0': 'x', '1': 'y', '2': 'z'}

        for kk, vv in axes.iteritems():

            prj = yt.ProjectionPlot(ds,
                                    int(kk),
                                    field,
                                    center='c')

            tc1 = [c for c in master_clump]
            prj.annotate_clumps(tc1)    # ok, seems like this works.., but we can't force it to plot other colors, okay, maybe the problem was the syntax.... it only takes
                            # plot_args={
                            # 'colors': [col_f(ff)]
                            # }

            # prj.set_unit(field=("h2density"), new_unit=..., equivalency=...)

            # print type(leaf_clumps)
            # print "***"
            # print type(tc1)
            # print "***"
            # import pdb; pdb.set_trace()
            prj.zoom(3)

            if saveplot:
                prj.save(fold_out + 'clumps1_' + '{0:.2f}'.format(n_cut) + \
                        '_' + str(int(step)) + '-' + str(int(N_cell_min)) + \
                        '_' + vv + 'axis.png')
            else:
                prj.show()

    if plot:
        plotclumps(ds, leaf_clumps, master_clump, saveplot=saveplot, fold_out=fold_out)

    return master_clump, leaf_clumps


if __name__ == '__main__':
    # -------------- decide on n_cut --------------
    # may be useful to plot the 3D, see test_rendering.py


    # Or see fig 6 Pallottini 2017 for n_H2 cut as starting point, lower right panel, based on Minkowsky function (Euler characteristic).

    # in units of nH2/cc
    n_cut_1 = 10**0.5
    n_cut_2 = 10**-1.5


    # -------------- run clump finder -------------

    master10, leaf10 = ytclumpfind_H2(ds, dd, ("h2density"),
                                      n_cut=0.1,
                                      step=5.0,
                                      N_cell_min=10,
                                      plot=True,
                                      saveplot=True,
                                      fold_out=fold_out)

    # write_clump_index(master10, 0, "master10_clump_hierarchy_H2.txt")
    # write_clumps(master10, 0,  "master10_clumps_H2.txt")

    # import os
    # os.system('cat *_clump_hierarchy_H2.txt')
    # os.system('cat *_clumps_H2.txt')


    # for ind in range(len(leaf10)):
    #     print(leaf10[ind]["h2density"])
    #     print(leaf10[ind].quantities.total_mass())
    #     print(leaf10[ind].quantities.center_of_mass())


    # master20, leaf20 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_1,
    #                                   step=20,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master20, 0, "master20_clump_hierarchy_H2.txt")
    # write_clumps(master20, 0,  "master20_clumps_H2.txt")

    # master30, leaf30 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_1,
    #                                   step=30,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master30, 0, "master30_clump_hierarchy_H2.txt")
    # write_clumps(master30, 0,  "master30_clumps_H2.txt")

    # master70, leaf70 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_1,
    #                                   step=70,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master70, 0, "master70_clump_hierarchy_H2.txt")
    # write_clumps(master70, 0,  "master70_clumps_H2.txt")

    # master100, leaf100 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_1,
    #                                     step=100,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # write_clump_index(master100, 0, "master100_clump_hierarchy_H2.txt")
    # write_clumps(master100, 0,  "master100_clumps_H2.txt")


    # master200, leaf200 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_1,
    #                                     step=200,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # write_clump_index(master200, 0, "master200_clump_hierarchy_H2.txt")
    # write_clumps(master200, 0,  "master200_clumps_H2.txt")


    # # --- repeat for n_cut_2 ---
    # master10, leaf10 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=10,
    #                                   N_cell_min=20,
    #                                   plot=True,
    #                                   saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master10, 0, "master10_clump_hierarchy_H2_cut2.txt")
    # write_clumps(master10, 0,  "master10_clumps_H2_cut2.txt")

    # master20, leaf20 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=20,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master20, 0, "master20_clump_hierarchy_H2_cut2.txt")
    # write_clumps(master20, 0,  "master20_clumps_H2_cut2.txt")

    # master30, leaf30 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=30,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master30, 0, "master30_clump_hierarchy_H2_cut2.txt")
    # write_clumps(master30, 0,  "master30_clumps_H2_cut2.txt")

    # master70, leaf70 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=70,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # try:
    #     write_clump_index(master70, 0, "master70_clump_hierarchy_H2_cut2.txt")
    #     write_clumps(master70, 0,  "master70_clumps_H2_cut2.txt")
    # except:
    #     pass

    # master100, leaf100 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_2,
    #                                     step=100,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # try:
    #     write_clump_index(master100, 0, "master100_clump_hierarchy_H2_cut2.txt")
    #     write_clumps(master100, 0,  "master100_clumps_H2_cut2.txt")
    # except:
    #     pass

    # master200, leaf200 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_2,
    #                                     step=200,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # try:
    #     write_clump_index(master200, 0, "master200_clump_hierarchy_H2_cut2.txt")
    #     write_clumps(master200, 0,  "master200_clumps_H2_cut2.txt")
    # except:
    #     pass

    # # find out leaf_clumps attributes to retreive physical properties
    # aa = leaf10[0]
    # print aa.field
    # print aa.data.fcoords    # grids position
    # print aa.data.icoords    # grids element index
    # print aa.data.icoords[0]
    # print aa.data['h2density'][0]
    # _density = f["rho"].value
    # _h2 = f["H2"].value
    # ii = aa.data.icoords[:, 0]
    # jj = aa.data.icoords[:, 1]
    # kk = aa.data.icoords[:, 2]
    # print (_density[ii, jj, kk] * factor) * _h2[ii, jj, kk]
    # assert round(aa.data['h2density'][0]) == round((_density[ii, jj, kk] * factor * _h2[ii, jj, kk])[0])

