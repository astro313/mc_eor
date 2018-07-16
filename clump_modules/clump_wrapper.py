"""
NOTE
----
The yt clump finder was initially described in http://adsabs.harvard.edu/abs/2009ApJ...691..441S , but it's changed since then. What it does now is decompose into non-overlapping tiles (stored in a kd-tree), identify contours within a tile, and then connect them across tiles. It does this with an upper and lower bound on a field value, and looking for topologically connected sets.

With yt, connected sets can be identified. This enables the location and analysis of a hierarchy of clumps; in star formation, for instance, this allows the construction of diagrams describing the density at which fragmentation occurs.

- extract_connected_sets() differ slightly from Clump()
    - extract_connected_sets(): get topologically connected sets. These sets are identified by examining cells between two threshold values and connecting them.
    - Clump(): uses a contouring algorithm to identified topologically disconnected structures within a dataset. This works by first creating a single contour over the full range of the contouring field, then continually increasing the lower value of the contour until it reaches the maximum value of the field. As disconnected structures are identified as separate contours, the routine continues recursively through each object, creating a hierarchy of clumps. Individual clumps can be kept or removed from the hierarchy based on the result of user-specified functions, such as checking for gravitational boundedness.

"""

import yt
from yt.analysis_modules.level_sets.api import Clump, find_clumps, get_lowest_clumps, write_clump_index, write_clumps

def ytclumpfind_H2(ds, dd, field, n_cut, c_max=None, step=3, N_cell_min=10, save=False, plot=False, saveplot=False, fold_out='./'):
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
        assert saveplot is True

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

    if plot:
        plotclumps(ds=ds, leaf_clumps = leaf_clumps,field=field, master_clump = master_clump, saveplot=saveplot, fold_out=fold_out
            , n_cut = n_cut,N_cell_min = N_cell_min,step = step
            )

    return master_clump, leaf_clumps

def plotclumps(ds, leaf_clumps, master_clump, field='h2density', saveplot=False, fold_out='',
    n_cut = 0 ,N_cell_min = 0,step= 0 ):
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



