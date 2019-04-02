#!/mnt/home/daisyleung/Downloads/yt-conda/bin/python

"""

Get to the buttom of the buggy yt clump finder...

Loop through all snapshots, plot h2density histogram, pick cuts, then find one clump at a time (vs. plot all clumps given a list of n_cuts, but there's a bug in yt that doesn't do this properly....)

NOTE
----
THE BUG (yt will reset itself to look from region > parent region...something like below):

min/max value for finding contours:  1.4677992676220695 92.60333032919448 g/cm**3
Number of children in this clump:  4
ncut =  1.4677992676220695
  N clumps 1
    id      0
    len     11239424
    max f*d 92.60333032919448 g/cm**3
    min f*d 8.966110228571374e-08 g/cm**3
Dump to
   test_brute/ss_28_ncut_1.47.png

-- Bug persist even if we do cut one-by-one. We modified yt clump_handling.py to fix this..


last mod: 23 July 2018


"""

import sys
assert 'yt' in sys.path[1]     # instead of using the anaconda version

import yt
import os
from yt.funcs import mylog
# This sets the log level to "ERROR", http://yt-project.org/doc/faq/index.html
mylog.setLevel(0)
import numpy as np

from clump_modules.clump_wrapper import ytclumpfind_H2, get_phyprop_of_leaf
from io_modules.manipulate_fetch_gal_fields import import_fetch_gal, import_fetch_stars, prepare_unigrid, check_hist_h2, check_power

from plot_modules.plot_cloud_prop import setup_plot
setup_plot()


outdir = 'test_brute/'

if not os.path.isdir(outdir):
    os.mkdir(outdir)

field_select = "h2density"

th_list = 10**np.linspace(-0.5, 1.5, 10)
#th_list = [31.62]

n_cell_min = 10
largeNum = 1.e+42   # to plot only one contour in a hacky way


istart = 16
iend = 29


isnap = str(istart)
test =  not os.path.isfile('snapshot'+isnap+'_center_stars_resampled.h5')
read_proper_unit = True

if read_proper_unit:

    from io_modules.load_misc import get_camera_from_file

    try:
      here = os.path.dirname(os.path.abspath(__file__))+'/'
    except NameError:
      # for ipython copy and paste compatibility
      here = '/mnt/home/daisyleung/mc_eor/'

    folder         = here+'precomputed_data/'
    f_camera       = folder + 'camera_settings.log'
    cameraDat = get_camera_from_file(f_camera)
else:
    pass

# I have some problem in printing the file (non unicode character somewhere)
debug = True


# because of the stupid yt bug, we will loop through the cuts and run
# clumpfinder one level at a time...

for snapshotnum in range(istart, iend):

    if read_proper_unit:
        regionsize_kpc = cameraDat[str(snapshotnum)]['size_kpc']
    else:
        regionsize_kpc = None

    # read data
    data = import_fetch_gal(isnap=snapshotnum)
    if not test:
      starData = import_fetch_stars(isnap=snapshotnum, verbose=False)

    ds, dd = prepare_unigrid(data=data,
                             add_unit=read_proper_unit,
                             regionsize_kpc=regionsize_kpc,
                             debug=False)

    check_hist_h2(data, th_list, ss=snapshotnum, outdir=outdir)
    # import pdb; pdb.set_trace()

    #check_power(data=data, size_kpc = 7., isnap=snapshotnum, outdir = outdir)

    # loop over cuts
    for incut in th_list:

        f_out = outdir + "ss_" + str(snapshotnum) + \
          "_ncut_" + "{:.2f}.pdf".format(incut)

        prj = yt.ProjectionPlot(ds, 0, field_select,
                                center='c', weight_field='h2density')

        #prj.set_zlim('h2density', 1.e-3, 1e-1)
        _, leaf_clumps = ytclumpfind_H2(ds, dd, field_select, incut,
                                        c_max=None, step=1e+6,
                                        N_cell_min=n_cell_min, save=False,
                                        plot=False, saveplot=None, fold_out='./')

        from plot_modules.dual_view_plotter import plot_face_edge

        if not debug:
          # plot face on and edgeon view of the galaxy and overplot clumps
          __ppj, __paxes = plot_face_edge(isnap=snapshotnum,
                                              data=data, ds=ds, dd=dd,
                                              leaf_clumps=leaf_clumps,
                                              incut=incut,
                                              selected_field='h2density',
                                              f_out=f_out.replace("ss", "dual"))
        # stupid yt does not check N_cell_min criterion if no children survived!
        if len(leaf_clumps) == 1:
            if (len(leaf_clumps[0]["h2density"]) < n_cell_min) or (len(leaf_clumps[0]["h2density"] > 100**3)):
                print("Removing the very last clump after checking N_cell_min criterion. Found no clumps...")
                del leaf_clumps
                continue
            if len(leaf_clumps[0]["h2density"]) == dd.shape[0]:
                print("Found no clumps and yt is being stupid..")
                del leaf_clumps
                continue

        id_sorted = range(len(leaf_clumps))
        # leaf_clumps to sort by np.sum(leaf_clumps[id_sorted]["density"])
        id_sorted = sorted(range(len(leaf_clumps)),
                           key=lambda x: np.sum(leaf_clumps[x]["density"]))

        prj.annotate_contour(field="h2density", ncont=1, factor=1,
                             clim=(incut, largeNum))  # to deal w/ stupid yt annotate_clump() bug

        # save clump positions
        outname = 'ss' + str(snapshotnum) + '_' + field_select + \
                  '_clumppos_ncut_' + '{:.2f}'.format(incut) + \
                  '_Ncellmin_' + str(int(n_cell_min)) + '.txt'
        X = []
        Y = []

        nn = int(len(dd['x'])**(1./3.)) + 1
        x_mat = np.copy(dd['x'].convert_to_units('kpc').value).reshape((nn,nn,nn))
        y_mat = np.copy(dd['y'].convert_to_units('kpc').value).reshape((nn,nn,nn))
        z_mat = np.copy(dd['z'].convert_to_units('kpc').value).reshape((nn,nn,nn))
        wg_mat= np.copy(dd['density'].convert_to_units('g/cm**3').value).reshape((nn,nn,nn))

        head  = 'Coordinates of clumps\nClumpID x, y, z [kpc]'
        head2 = 'Weighted coordinates of clumps\nClumpID x, y, z [kpc]'
        for ileaf in id_sorted:
            # unweighted center
            _fc = np.mean(leaf_clumps[ileaf].data.fcoords[:], axis=0)
            # weighted center
            ids = leaf_clumps[ileaf].data.icoords
            ix,iy,iz  = ids[:,0], ids[:,1], ids[:,2]
            _fc_wg    = np.array([0.,0.,.0])
            _fc_wg[0] = np.sum(wg_mat[ix,iy,iz]*x_mat[ix,iy,iz])/np.sum(wg_mat[ix,iy,iz])
            _fc_wg[1] = np.sum(wg_mat[ix,iy,iz]*y_mat[ix,iy,iz])/np.sum(wg_mat[ix,iy,iz])
            _fc_wg[2] = np.sum(wg_mat[ix,iy,iz]*z_mat[ix,iy,iz])/np.sum(wg_mat[ix,iy,iz])

            if debug:
              print ileaf
              print '  c_wg',_fc_wg
              print '  c   ',_fc
              print '  wg r',np.max(wg_mat[ix,iy,iz]),np.min(wg_mat[ix,iy,iz])

            X.append([ileaf, _fc[0], _fc[1], _fc[2]])
            Y.append([ileaf, _fc_wg[0], _fc_wg[1], _fc_wg[2]])
        print 'save to'
        print '  ',outname
        np.savetxt(outname, X, fmt='%d %.5f %.5f %.5f', header=head)
        np.savetxt(outname.replace('_clumppos_','_wgclumppos_'), Y, fmt='%d %.5f %.5f %.5f', header=head2)

        for ileaf in id_sorted:
            _fc = np.mean(leaf_clumps[ileaf].data.fcoords[:], axis=0)
            prj.annotate_marker(_fc,
                                coord_system='data',
                                plot_args={'color': 'red', 's': 500})
            prj.annotate_text(_fc,
                              ileaf,
                              coord_system='data',
                              text_args={'color': 'black', 'size': 32},
                              inset_box_args={'boxstyle': 'square',
                                              'facecolor': 'white',
                                              'linewidth': 2.0,
                                              'edgecolor': 'white',
                                              'alpha': 0.})

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
        prj = None

        import cPickle as pickle
        print '  compute properties'
        # to retreive physical properties of leaf
        leaf_fields = {}
        for n_leaf in range(len(leaf_clumps)):
            leaf_fields[str(n_leaf)] = get_phyprop_of_leaf(leaf_clumps[n_leaf],
                                                           data['density'],
                                                           data['H2'] *
                                                           data['density'],
                                                           data['P'], data[
                                                               'P_nt'],
                                                           data['Z'],
                                                           data['velx'],
                                                           data['vely'],
                                                           data['velz'],
                                                           starPartDict=starData,
                                                           plothist=False)
        print "saved leaf fields: ", leaf_fields['0'].keys()
        suboutdir = "leaf_fields_" + str(snapshotnum) + "/"

        if not os.path.isdir(outdir + suboutdir):
            os.mkdir(outdir + suboutdir)

        pickle.dump(leaf_fields, open(outdir + suboutdir +
                                      '{0:.2f}'.format(incut) + '_' + str(n_cell_min) + "_fields.p", "wb"))

    # ------------
