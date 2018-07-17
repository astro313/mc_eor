
"""

Example usage:
    python snapshot28_ytclump_fields.py 28 --savepickle --plot --saveplot test_png/


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
import h5py

from yt.analysis_modules.level_sets.api import Clump

def col_f(ii, cm=None):
    """
    get a different color for each indx, ii
    """
    if cm is None:
        cm = plt.get_cmap('gist_heat')
    return cm(ii)

from clump_modules.clump_wrapper import ytclumpfind_H2

def get_phyprop_of_leaf(subleaf, density, H2density, Pressure, P_nt, metallicity, velx, vely, velz, plothist=False):
    """

    Parameters
    ----------
    subleaf: a leaf object (i.e., w/o children)

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

    # index corresponding to the clumps (w/ no children)
    ii = subleaf.data.icoords[:, 0]
    jj = subleaf.data.icoords[:, 1]
    kk = subleaf.data.icoords[:, 2]

    # print(subleaf.quantities.total_mass())
    # print(subleaf.quantities.center_of_mass())

    for fff in _leaf_fields.iterkeys():
        _leaf_fields[fff] = eval(fff)[ii, jj, kk]

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

    return _leaf_fields


if __name__ == '__main__':

    import argparse
    import textwrap
    import os

    desc = """ find clumps/molecular complexes from unigrid cube. """

    epi = textwrap.dedent('''
        The clump-find is done using yt, which decomposes into non-overlapping tiles (stored in a kd-tree), identify contours within a tile, and then connect them across tiles. It does this with an upper and lower bound on a field value, and looking for topologically connected sets.

        The fit parameters are (in this order):

            snapshot_num: snapshot file to load in, integer.

            not_convert_unit: if true, convert NOT from code unit to more commonly used units, depending on how fetch_gal_fields.py is implemented.

            n_cut: threshold to look for clumps, in units of nH2/cc

            step: multiplicative interval between subsequent contours

            N_cell_min: min. number of cell s.t. clump identified is not spurious

            save: if true, save the clump tree as a reloadable dataset (using yt func)

            savepickle: if true, save the fields of all leafs stored in dict into a pickled file

            plot: if true, will plot leaf clump

            saveplot: if true, will save figure instead of showing it

            fold_out: directory to save figures

        ''')

    fmter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=fmter)

    # -------------- decide on n_cut --------------
    # fig 6 Pallottini 2017 for n_H2 cut as starting point, lower right panel,
    # based on Minkowsky function (Euler characteristic).

    # in units of nH2/cc
    # n_cut_1 = 10**0.5
    n_cut_2 = 10**-1.5

    # -------------- parse arguments ----------------

    parser.add_argument('snapshot_num', action="store", type=int,
                        help="snapshot number to load in (no default).")

    parser.add_argument('--not_convert_unit',
                        action="store_true",
                        default=False,
                        help="do NOT convert from code units to more commonly used units, depending on how fetch_gal_fields.py is implemented.")

    parser.add_argument('-nc', '--ncut', action="store", type=float,
                        default=n_cut_2,
                        help="threshold to look for clumps, in units of nH2 [1/cc]")

    parser.add_argument('-s', '--step', action="store", type=int,
                        default=5,
                        help="multiplicative interval between subsequent contours")

    parser.add_argument('-nm', '--Nmin', action="store", type=int,
                        default=3,
                        help="min. number of cell s.t. clump identified is not spurious")

    parser.add_argument('--save', action="store_true", default=False,
                        help="save the clump tree as a reloadable dataset (using yt func)")

    parser.add_argument('--savepickle', action="store_true", default=False,
                        help="save the fields of all leafs stored in dict into a pickled file")

    parser.add_argument('--plot', action="store_true", default=False,
                        help="plot leaf clump on projection plots")

    parser.add_argument('--saveplot', action="store_true", default=False,
                        help="save figure instead of showing it")

    parser.add_argument('fold_out', action='store',
                        default='test_png/',
                        help="directory to save plot files ending with /")

    parser.add_argument('--verbose', action="store_true", default=False,
                        help="verbose option")

    args = parser.parse_args()

    # ---------------------------------------------------------------


    from io_modules.manipulate_fetch_gal_fields import import_fetch_gal, prepare_unigrid

    data = import_fetch_gal(isnap = args.snapshot_num, verbose = args.verbose , convert = (not args.not_convert_unit))

    ds,dd = prepare_unigrid(data= data)


    # -------------- run clump finder -------------
    # master5, leaf5 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                 n_cut=n_cut_1,
    #                                 step=5,
    #                                 N_cell_min=3,
    #                                 plot=True,
    #                                 saveplot=True,
    #                                 fold_out=fold_out)

    assert isinstance(args.fold_out, str)
    if not os.path.isdir(args.fold_out):
        os.mkdir(args.fold_out)

    if(args.verbose):
        print 'start clumpfinder'
    # --- repeat for n_cut_2 ---
    master5, leaf5 = ytclumpfind_H2(ds, dd, ("h2density"),
                                    n_cut=args.ncut,
                                    step=args.step,
                                    N_cell_min=args.Nmin,
                                    plot=args.plot,
                                    saveplot=args.saveplot,
                                    fold_out=args.fold_out)
    if(args.verbose):
        print '  complete'

    if(args.verbose):
        print '  compute properties'
    # to retreive physical properties of leaf
    leaf_fields = {}
    for n_leaf in range(len(leaf5)):
        leaf_fields[str(n_leaf)] = get_phyprop_of_leaf(leaf5[n_leaf],
                                                       data['density'],
                                                       data['H2'] * data['density'],
                                                       data['P'], data['P_nt'],
                                                       data['Z'],
                                                       data['velx'], data['vely'], data['velz'],
                                                       plothist=False)
    print "saved leaf fields: ", leaf_fields['0'].keys()

    if args.savepickle:
        import cPickle as pickle
        # import pickle

        outdir = "leaf_fields_" + str(args.snapshot_num) + "/"

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        pickle.dump(leaf_fields, open(outdir + '{0:.2f}'.format(args.ncut) + '_' + str(
            args.step) + '_' + str(args.Nmin) + "_fields.p", "wb"))


# -------------
