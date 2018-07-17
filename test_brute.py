"""

Get to the buttom of the buggy yt clump finder...

Loop through all snapshots, plot h2density histogram, pick cuts, then find one clump at a time (vs. plot all clumps given a list of n_cuts, but there's a bug in yt that doesn't do this properly....)

NOTE
----
THE BUG (yt will reset itself to look from region > parent region...something like below):
Number of children in this clump:  8
ncut =  14.67799267622069
  N clumps 1
    id      0
    len     16777216
    max f*d 48.6998896743599 g/cm**3
    min f*d 6.880567814353917e-08 g/cm**3


last mod: 16 July 2018


"""

import matplotlib.pyplot as plt
import yt, os
from yt.funcs import mylog
# This sets the log level to "ERROR", http://yt-project.org/doc/faq/index.html
mylog.setLevel(40)
import numpy as np

from clump_modules.clump_wrapper import ytclumpfind_H2
from io_modules.manipulate_fetch_gal_fields import import_fetch_gal, prepare_unigrid, check_hist_h2


def col_f(ii, cm=None):
    if cm is None:
        cm = plt.get_cmap('gist_heat')
    return cm(ii)


outdir = 'test_brute/'

if not os.path.isdir(outdir):
    os.mkdir(outdir)

field_select = "h2density"

th_list = 10**np.linspace(-0.5, 1.5, 7)
# th_list = [7.]
# th_list = [10**-0.5]
test = True

n_cell_min = 10

# because of the stupid yt bug, we will loop through the cuts and run clumpfinder one level at a time...
for incut in th_list:

    # loop through all snapshots
    for snapshotnum in range(16, 29):
        data = import_fetch_gal(isnap=snapshotnum)

        check_hist_h2(data, incut, ss=snapshotnum, outdir=outdir)
        f_out = outdir + "ss_" + str(snapshotnum) + "_ncut_" + "{:.2f}.png".format(incut)

        ds, dd = prepare_unigrid(data=data)
        prj = yt.ProjectionPlot(ds, 0, field_select,
                                center='c', weight_field='density')

        #prj.set_zlim('h2density', 1.e-3, 1e-1)
        #th_list = np.linspace(0.01, .1, 3)
        # th_list = [0.01]

        _, leaf_clumps = ytclumpfind_H2(ds, dd, field_select, incut,
                                            c_max=None, step=1e+6,
                                            N_cell_min=n_cell_min, save=False,
                                            plot=False, saveplot=None, fold_out='./')

        id_sorted = range(len(leaf_clumps))
        # leaf_clumps to sort by np.sum(leaf_clumps[id_sorted]["density"])
        id_sorted = sorted(range(len(leaf_clumps)),
                               key=lambda x: np.sum(leaf_clumps[x]["density"]))

        for ileaf in id_sorted:
            _fc = np.mean(leaf_clumps[ileaf].data.fcoords[:], axis=0)

            if(0):
                ff = float(iii) / len(th_list)
            else:
                ff = float(ileaf) / len(leaf_clumps)

            prj.annotate_clumps([leaf_clumps[ileaf]], plot_args={
                                    'colors': [col_f(ff)],
                                    'alpha': 0.8,
                                    'zorder': ileaf
                                    })

            prj.annotate_marker(_fc,
                                    coord_system='data',
                                    plot_args={'color': 'red', 's': 500})
            prj.annotate_text(_fc,
                                  ileaf,
                                  coord_system='data',
                                  text_args={'color': 'black', 'size': 8},
                                  inset_box_args={'boxstyle': 'square',
                                                  'facecolor': 'white',
                                                  'linewidth': 2.0,
                                                  'edgecolor': 'white',
                                                  'alpha': 0.35})

        print 'ncut = ', incut
        print '  N clumps', len(leaf_clumps)
        for i in xrange(len(leaf_clumps)):
            print '    id     ', i
            print '    len    ', len(leaf_clumps[i]["h2density"])
            print '    max f*d', np.max(leaf_clumps[i]["h2density"])
            print '    min f*d', np.min(leaf_clumps[i]["h2density"])

        prj.zoom(2)
        prj.save(f_out)
        print 'Dump to'
        print '  ', f_out
        print '\n'

        prj = None


        from snapshot28_ytclump_fields import get_phyprop_of_leaf
        import cPickle as pickle
        print '  compute properties'
        # to retreive physical properties of leaf
        leaf_fields = {}
        for n_leaf in range(len(leaf_clumps)):
            leaf_fields[str(n_leaf)] = get_phyprop_of_leaf(leaf_clumps[n_leaf],
                                                       data['density'],
                                                       data['H2'] * data['density'],
                                                       data['P'], data['P_nt'],
                                                       data['Z'],
                                                       data['velx'],
                                                       data['vely'],
                                                       data['velz'],
                                                       plothist=False)
        print "saved leaf fields: ", leaf_fields['0'].keys()
        suboutdir = "leaf_fields_" + str(snapshotnum) + "/"

        pickle.dump(leaf_fields, open(outdir + suboutdir + '{0:.2f}'.format(ncut) + '_' + str(n_cell_min) + "_fields.p", "wb"))


    # ------------
