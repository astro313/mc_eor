'''

plot physical properties of clouds across all snapshots, identified using same way (ncut, steps, etc).

last mod: 24 July 2018


'''

import os
from plot_modules.plot_cloud_prop import setup_plot
plt = setup_plot()


def plotting_procedure(minss, maxss, fname):
    from plot_modules.plot_cloud_prop import unpack_xy, plot_stuff, plot_stuff_3dim
    from io_modules.leaf_pickle import load_in_pickle_SS

    ss = load_in_pickle_SS('/mnt/home/daisyleung/mc_eor/test_brute/leaf_fields_', fname, minss, maxss)
    to_plot = unpack_xy(ss)

    outleafdir = '/mnt/home/daisyleung/mc_eor/test_brute/leaf_fields_ss' + \
                 str(minss) + '-' + str(maxss) + '/'
    leafdir_out = outleafdir + fname[:fname.find('.p')] + '/'
    if not os.path.isdir(leafdir_out):
        os.system('mkdir -p ' + leafdir_out)

    plot_stuff("cloud mass", "mass over jeans mass", leglabel="sfr: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "jeans mass",  leglabel="sfr: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "alpha vir", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("size pc", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("Mach", "SFR young", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("stellar to gas mass", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
#     plot_stuff("gas sd", "sfr sd", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("gas sd", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("R2 pc2", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff_3dim("tff Myr", "size pc", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff('gas sd cgs', 'sigmaSq over size', leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)

    from plot_modules.plot_cloud_prop import get_masses_all_clouds, massFuncUnbinnedCDF, massFuncPDF, massFuncDifferential
    allmass = get_masses_all_clouds(ss)
    massFuncUnbinnedCDF(allmass, outdir=leafdir_out)
    massFuncPDF(allmass, outdir=leafdir_out)
    massFuncDifferential(allmass, outdir=leafdir_out)


# ---

if __name__ == '__main__':
    plotting_procedure(16, 28, fname="0.32_10_fields.p")
    plotting_procedure(16, 28, fname="0.53_10_fields.p")
    plotting_procedure(16, 28, fname="0.88_10_fields.p")
    plotting_procedure(16, 28, fname="1.47_10_fields.p")
    plotting_procedure(16, 28, fname="2.45_10_fields.p")
    plotting_procedure(16, 28, fname="4.08_10_fields.p")
    plotting_procedure(16, 28, fname="6.81_10_fields.p")
    plotting_procedure(16, 28, fname="11.36_10_fields.p")
    plotting_procedure(16, 28, fname="18.96_10_fields.p")
