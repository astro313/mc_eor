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
    # n_cut = 1.0    # to make sure whatever comes out after multiplicative by
    # step won't be too small

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
        plotclumps(ds=ds, leaf_clumps=leaf_clumps, field=field, master_clump=master_clump, saveplot=saveplot, fold_out=fold_out, n_cut=n_cut, N_cell_min=N_cell_min, step=step
                   )

    return master_clump, leaf_clumps


def plotclumps(ds, leaf_clumps, master_clump, field='h2density', saveplot=False, fold_out='',
               n_cut=0, N_cell_min=0, step=0):
    """ overplot the clumps found (specifically the leaf_clumps) along 3 images, each created by projecting onto x-, y-, and z-axis. """

    axes = {'0': 'x', '1': 'y', '2': 'z'}

    for kk, vv in axes.iteritems():

        prj = yt.ProjectionPlot(ds,
                                int(kk),
                                field,
                                center='c')

        tc1 = [c for c in master_clump]
        # ok, seems like this works.., but we can't force it to plot other
        # colors, okay, maybe the problem was the syntax.... it only takes
        # plot_args={
        # 'colors': [col_f(ff)]
        # }
        prj.annotate_clumps(tc1)

        # print type(leaf_clumps)
        # print "***"
        # print type(tc1)
        # print "***"
        # import pdb; pdb.set_trace()
        prj.zoom(3)

        if saveplot:
            prj.save(fold_out + 'clumps1_' + '{0:.2f}'.format(n_cut) +
                     '_' + str(int(step)) + '-' + str(int(N_cell_min)) +
                     '_' + vv + 'axis.png')
        else:
            prj.show()


def get_phyprop_of_leaf(subleaf, density, H2density, Pressure, P_nt, metallicity, velx, vely, velz, starPartDict=None, plothist=False):
    """

    Parameters
    ----------
    subleaf: a "leaf" object

    density: array
        density of uniformly sampled subregion (probably converted into desire units at this point.)

    H2density: array
        H2 density of uniformly sampled subregion (probably converted into desire units at this point.)

    Pressure: array
        Thermal Pressure of uniformly sampled subregion (probably converted into desire units at this point.)

    P_nt: array
        P_nt of uniformly sampled subregion

    metallicity: array
        Z of uniformly sampled subregion, divide by 0.02 to get to solar metallicity

    velx: array
        vel-x of uniformly sampled subregion

    vely: array
        vel-y of uniformly sampled subregion

    velz: array
        vel-z of uniformly sampled subregion

    starPartDict: dict
        contains mass and epochs fields of stars

    plothist: bool
        whether or not to plot the histogram for each field corresponding to the clump

    Returns
    -------
    _leaf_fields: dict
        bunch of fields

    """

    _leaf_fields = {'density': None,
                    'H2density': None,
                    'Pressure': None,
                    'P_nt': None,
                    'metallicity': None,
                    'velx': None,
                    'vely': None,
                    'velz': None}

    if starPartDict is not None:
        for kkk in starPartDict.iterkeys():
            _leaf_fields[kkk] = None

    # index corresponding to the clumps
    ii = subleaf.data.icoords[:, 0]
    jj = subleaf.data.icoords[:, 1]
    kk = subleaf.data.icoords[:, 2]

    # print(subleaf.quantities.total_mass())
    # print(subleaf.quantities.center_of_mass())

    for fff in _leaf_fields.iterkeys():
        if not fff in ['mass', 'epoch']:
            _leaf_fields[fff] = eval(fff)[ii, jj, kk]
        else:
            _leaf_fields[fff] = starPartDict[kkk][ii, jj, kk]

    #_h2_subleaf = subleaf.data['h2density']
    #h2_subleaf = _leaf_fields['H2density']
    #assert round(_h2_subleaf[0]) == round(h2_subleaf[0])
    #del h2_subleaf
    #del _h2_subleaf

    if plothist:
        # hard code for now
        plt.hist(density_aa)
        plt.title('density [1/cm^3]')
        plt.show()
        plt.hist(H2_aa)
        plt.title('H2 density [1/cm^3]')
        plt.show()
        plt.hist(np.log10(Pressure_aa))
        plt.title('Thermal Pressure [K cm-3]')
        plt.show()
        plt.hist(np.log10(P_nt_aa))
        plt.title('Non-thermal Pressure [K cm-3]')
        plt.show()
        plt.hist(metallicity_aa / 0.02)
        plt.title('Metallicity in Solar Z units')
        plt.show()
        plt.hist(velx)
        plt.title('vel-x [km/s]')
        plt.show()
        plt.hist(vely)
        plt.title('vel-y [km/s]')
        plt.show()
        plt.hist(velz)
        plt.title('vel-z [km/s]')
        plt.show()

        if starPartDict is not None:
            plt.hist(np.log10(mass[mass > 1.e-23]))
            plt.title('mass')
            plt.show()

    return _leaf_fields
