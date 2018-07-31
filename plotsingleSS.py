'''

plot physical properties of clouds.

Single snapshot, different ncuts (from test_brute.py)


last mod: 30 July 2018

'''


import os
from plot_modules.plot_cloud_prop import setup_plot
plt = setup_plot()

def plotting_procedure(snapshot_num):

    try:
      here = os.path.dirname(os.path.abspath(__file__))+'/'
    except NameError:
      # for ipython copy and paste compatibility
      here = '/mnt/home/daisyleung/mc_eor/'

    leafdir_out = here+"test_brute/leaf_fields_" + str(snapshot_num) + '/'

    if not os.path.isdir(leafdir_out):
        os.mkdir(leafdir_out)

    from io_modules.leaf_pickle import load_in_pickled_leaf_singleSS
    from plot_modules.plot_cloud_prop import unpack_xy, plot_stuff, plot_stuff_3dim

    ss = load_in_pickled_leaf_singleSS(leafdir_out, snapshot_num)
    to_plot = unpack_xy(ss)

    plot_stuff("cloud mass", "mass over jeans mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff("cloud mass", "jeans mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, cm='gist_rainbow')
    plot_stuff("cloud mass", "alpha vir", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff("size pc", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff("Mach", "SFR young", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff("stellar to gas mass", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
#     plot_stuff("gas sd", "sfr sd", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff("gas sd", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff("R2 pc2", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff_3dim("tff Myr", "size pc", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
    plot_stuff('gas sd cgs', 'sigmaSq over size', leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)


    from plot_modules.plot_cloud_prop import get_masses_all_clouds, massFuncUnbinnedCDF, massFuncPDF, massFuncDifferential
    allmass = get_masses_all_clouds(ss)
    tag = 'ss'+str(snapshot_num) +'diffncuts'
    massFuncUnbinnedCDF(allmass, outdir=leafdir_out, tag=tag)
    massFuncPDF(allmass, outdir=leafdir_out, tag=tag)
    massFuncDifferential(allmass, outdir=leafdir_out, tag=tag)


# ---

if __name__ == '__main__':

    for isnap in range(16, 29):
        plotting_procedure(isnap)

