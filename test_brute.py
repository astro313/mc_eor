import h5py
import yt
from yt.funcs import mylog
mylog.setLevel(40) # This sets the log level to "ERROR", http://yt-project.org/doc/faq/index.html
import numpy as np
from test_refined_snapshot28_ytclump_fields import ytclumpfind_H2

field = ("h2density")

f = h5py.File("snapshot" + str(28) +
              "_center_fields0123456-15_resampled.h5", "r")
# careful sometimes i used "density" (e.g., resample.py), see
# resample_fields.py to make sure
density = f["rho"].value
H2 = f["H2"].value
Pressure = f["P"].value
P_nt = f["P_nt"].value
metallicity = f["Z"].value
velx = f["vel_x"].value
vely = f["vel_y"].value
velz = f["vel_z"].value


import pymses
from pymses.utils import constants as C

ro = pymses.RamsesOutput("output", 28)

from fetch_gal_fields import get_units
factor_density = get_units(ro=ro)['rho'][0]      # 1/cm^3 (not H/cm^3)
density *= factor_density
print density.max()

factor_vel = get_units(ro=ro)['vel'][0]
velx *= factor_vel
vely *= factor_vel
velz *= factor_vel
print velx.max(), vely.max(), velz.max()

factor_P = get_units(ro=ro)['P'][0]
Pressure *= factor_P
P_nt *= factor_P
print np.log10(Pressure.max()), np.log10(P_nt.max())

# th_list = [1.e-2, 0.1, 0.5]
# th_list = [0.4]
# th_list = np.linspace(0.4,0.1,3)
th_list = 10**np.linspace(-0.5, 1.5, 7)

data = dict(density=density, H2=H2,
            P=Pressure,
            P_nt=P_nt,
            Z=metallicity,
            velx=velx,
            vely=vely,
            velz=velz
            )

#data["density"] = np.ones_like(data["H2"])
data["density"][data["density"]<=0] = np.min(data["density"][data["density"]>0])
data["H2"][data["H2"] < 1.e-3] = 1.e-3

for var in ['density', 'H2']:
    print '  ', var, np.max(data[var]), np.min(data[var])

import matplotlib.pyplot as plt
aa = (data["density"] * data["H2"]).flatten()
print np.max(aa), np.min(aa), np.min(aa[aa > 0])
aa[aa <= 0] = np.min(aa[aa > 0])
aa = np.log10(aa)

print np.max(aa), np.min(aa), np.min(aa[aa > 0])

plt.close('all')
plt.figure()
plt.hist(aa, bins=100)
for ele in th_list:
    x = np.log10(ele)
    plt.plot([x, x], [1, 1.e+7], ls='--', color='k')
plt.yscale('log')
plt.savefig('hist_test.png')

#raise TypeError


def _h2density(field, data):
    try:
        return data["density"] * data["H2"]
    except:
        return data[("stream", "density")] * data[("stream", "H2")]

from yt.units import dimensions
# The global yt.add_field() function is for adding a field for every
# subsequent dataset that is loaded in a particular python session,
# whereas ds.add_field() (add_field()) will only add it to dataset ds.

ds = yt.load_uniform_grid(data, f["rho"].shape)
dd = ds.all_data()
ds.add_field(("stream", "h2density"), function=_h2density, units="g/cm**3")
assert (dd['H2'] * dd['density']).max() == dd['h2density'].max()

prj = yt.ProjectionPlot(ds, 0,
                        field,
                        center='c', weight_field='density')

#prj.set_zlim('h2density', 1.e-3, 1e-1)
#th_list = np.linspace(0.01, .1, 3)
# th_list = [0.01]


def col_f(ii, cm=None):
    if cm is None:
        cm = plt.get_cmap('gist_heat')
    return cm(ii)

th_list = th_list[1:]

for iii, nthres in enumerate(th_list):

    cmax = None


    _, leaf_clumps = ytclumpfind_H2(ds, dd, field, nthres,
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
prj.save('12345.png')

prj = None

# ------------
