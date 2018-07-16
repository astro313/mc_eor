
import matplotlib.pyplot as plt
import yt
from yt.funcs import mylog
mylog.setLevel(40) # This sets the log level to "ERROR", http://yt-project.org/doc/faq/index.html
import numpy as np

from test_refined_snapshot28_ytclump_fields import ytclumpfind_H2
from io_modules.manipulate_fetch_gal_fields import import_fetch_gal, prepare_unigrid,check_hist_h2

def col_f(ii, cm=None):
    if cm is None:
        cm = plt.get_cmap('gist_heat')
    return cm(ii)


f_out        = '12345.png'

field_select = "h2density"
th_list      = 10**np.linspace(-0.5, 1.5, 7)

test         = False




data = import_fetch_gal(isnap = 28)

if test:
  check_hist_h2(data=data)

ds,dd = prepare_unigrid(data= data)

prj = yt.ProjectionPlot(ds, 0,field_select,
                        center='c', weight_field='density')

#prj.set_zlim('h2density', 1.e-3, 1e-1)
#th_list = np.linspace(0.01, .1, 3)
# th_list = [0.01]

for iii, nthres in enumerate(th_list):

    cmax = None


    _, leaf_clumps = ytclumpfind_H2(ds, dd, field_select, nthres,
                                    c_max=cmax, step=1e+6,
                                    N_cell_min=10, save=False,
                                    plot=False, saveplot=None, fold_out='./')
    # prj.annotate_clumps(leaf_clumps)       # this is a bug.. instead, loop
    # through each leaf_clump an plot the contours..

    id_sorted = range(len(leaf_clumps))
    #leaf_clumps to sort by np.sum(leaf_clumps[id_sorted]["density"])

    id_sorted = sorted(range(len(leaf_clumps)), key= lambda x: np.sum(leaf_clumps[x]["density"]))

    for ileaf in id_sorted:
        _fc = np.mean(leaf_clumps[ileaf].data.fcoords[:], axis=0)

        if(1):
          ff = float(iii) / len(th_list)
        else:
          ff = float(ileaf) / len(leaf_clumps)

        prj.annotate_clumps([leaf_clumps[ileaf]], plot_args={
                            'colors': [col_f(ff)],
                            'alpha': 0.8,
                            'zorder': ileaf
                            })
        '''
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
        '''
    print 'ncut = ', nthres
    print '  N clumps', len(leaf_clumps)
    for i in xrange(len(leaf_clumps)):
        print '    id     ', i
        print '    len    ', len(leaf_clumps[i]["h2density"])
        print '    max f*d', np.max(leaf_clumps[i]["h2density"])
        print '    min f*d', np.min(leaf_clumps[i]["h2density"])

prj.zoom(2)
prj.save(f_out)
print 'Dump to'
print '  ',f_out

prj = None

# ------------
