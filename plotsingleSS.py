'''

plot physical properties of clouds.

Single snapshot, different ncuts (from brute.py)


main plots:
- alpha-vir_cloud-mass.png
- sigma-kms_size-pc.png
- sigma-kms_stellar-to-gas-mass.png
- Mach_cloud-mass.png
- weighted-density_size-pc.png
- sigma-kms-bulk_sigma-kms-NT.png


last mod: August 18 2018

'''

def load_pickleTOplot(snapshot_num, pattern=None):

    from io_modules.leaf_pickle import load_in_pickled_leaf_singleSS
    from plot_modules.plot_cloud_prop import unpack_xy
    import os

    try:
        here = os.path.dirname(os.path.abspath(__file__)) + '/'
    except NameError:
        # for ipython copy and paste compatibility
        here = '/mnt/home/daisyleung/mc_eor/'

    leafdir_out = here + "test_brute/leaf_fields_" + str(snapshot_num) + '/'

    if not os.path.isdir(leafdir_out):
        os.mkdir(leafdir_out)

    if pattern is None:
        ss = load_in_pickled_leaf_singleSS(leafdir_out, snapshot_num)
    else:
        ss = load_in_pickled_leaf_singleSS(leafdir_out, snapshot_num, pattern=pattern)       # , pattern='2.45_*10*.p')

    to_plot = unpack_xy(ss)
    return ss, to_plot, leafdir_out


def plotting_procedure(snapshot_num, legendFontSize=13):
    from plot_modules.plot_cloud_prop import plot_stuff, plot_stuff_3dim

    ss, to_plot, leafdir_out = load_pickleTOplot(snapshot_num)

    plot_stuff("stellar to gas mass", "cloud mass", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("stellar to gas mass", "alpha vir", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("stellar to gas mass", "Mach", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("gas sd per ff", "sfr sd", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff('gas sd cgs', 'sigmaSq over size', leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff('alpha vir', 'sigmaSq over size', leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("alpha vir", "SFR young", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("cloud mass", "mass over jeans mass",
               leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("cloud mass", "jeans mass", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("cloud mass", "alpha vir", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    # alpha_vir w/ different definition of sigma in it
    plot_stuff("cloud mass", "alpha vir", compareAlphaVir=True,
               leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("size pc", "sigma kms", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("Mach", "SFR young", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("Mach", "SFR old", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("stellar to gas mass", "sigma kms", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("gas sd", "sigma kms", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("R2 pc2", "cloud mass", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("size pc", "cloud mass", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=12)
    plot_stuff("cloud mass", "Mach", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff("size pc", "weighted density",  leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)
    plot_stuff('sigma kms NT', 'sigma kms bulk', leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)

    plot_stuff_3dim("tff Myr", "size pc", "cloud mass",
                    leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, legendFontSize=legendFontSize)

    # from plot_modules.plot_cloud_prop import get_masses_all_clouds, massFuncUnbinnedCDF, massFuncPDF, massFuncDifferential
    # allmass = get_masses_all_clouds(ss)
    # tag = 'ss' + str(snapshot_num) + 'diffncuts'
    # massFuncUnbinnedCDF(allmass, outdir=leafdir_out, tag=tag)
    # massFuncPDF(allmass, outdir=leafdir_out, tag=tag)
    # massFuncDifferential(allmass, outdir=leafdir_out, tag=tag)



# ---
if __name__ == '__main__':

    import matplotlib
    # matplotlib.use('agg')

    from plot_modules.plot_cloud_prop import setup_plot
    cm = setup_plot()

    for isnap in range(16, 29):
        plotting_procedure(isnap)

    # # min MC mass for highest n_cut (for paper)
    # for isnap in [16, 27]:
    #     _, _, _ = load_pickleTOplot(isnap, pattern='*18.96*10*p')

    # # plot Larson and local points only (for talk introduction slide)
    # _, to_plot, _ = load_pickleTOplot(16)
    # plot_stuff("size pc", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir='./')

    # for paper, plot alpha_vir - M_cl for ss21 (pre-SB) and ss22
    import matplotlib.pyplot as plt
    plt.close('all')

    from plotsingleSS import load_pickleTOplot
    from plot_modules.plot_cloud_prop import plot_alpha_vir_2ss, plot_Mach_massRatio_4ss, plot_alpha_vir_8ss
    ss1, to_plot1, leafdir_out1 = load_pickleTOplot(21)
    ss2, to_plot2, leafdir_out2 = load_pickleTOplot(22)
    fig, ax = plot_alpha_vir_2ss(to_plot1, to_plot2, ls='',
                              markersize=10,
                              marker='*',
                              tag='ss21-ss22',
                              t1='Pre-starburst Phase',
                              t2='Starburst Phase',
                              cbarLabelSize=20,
                              outdir='./',
                              legendFontSize=16,
                              saveFig=True)

    # 4x2 alpha_vir - cloud mass
    ss1, to_plot1, leafdir_out1 = load_pickleTOplot(16)
    ss2, to_plot2, leafdir_out2 = load_pickleTOplot(17)
    ss3, to_plot3, leafdir_out3 = load_pickleTOplot(18)
    ss4, to_plot4, leafdir_out4 = load_pickleTOplot(19)
    ss5, to_plot5, leafdir_out5 = load_pickleTOplot(20)
    ss6, to_plot6, leafdir_out6 = load_pickleTOplot(21)  # Pre-SB
    ss7, to_plot7, leafdir_out7 = load_pickleTOplot(22)  # SB
    ss8, to_plot8, leafdir_out8 = load_pickleTOplot(23)
    fig, ax = plot_alpha_vir_8ss(to_plot1, to_plot5, to_plot2, to_plot6,
                                     to_plot3, to_plot7, to_plot4, to_plot8,
                                     ls='',
                                     markersize=10,
                                     marker='*',
                                     tag='ss16-ss23',
                                     cbarLabelSize=20,
                                     outdir='./',
                                     showcbar=True,
                                     legendFontSize=16,
                                     saveFig=True)

    # Mach - stellar-to-gas mass ratio
    ss1, to_plot1, leafdir_out1 = load_pickleTOplot(21)
    ss2, to_plot2, leafdir_out2 = load_pickleTOplot(22)   # SB
    ss3, to_plot3, leafdir_out3 = load_pickleTOplot(23)
    ss4, to_plot4, leafdir_out4 = load_pickleTOplot(27)
    fig, ax = plot_Mach_massRatio_4ss(to_plot1, to_plot2, to_plot3, to_plot4,
                                      ls='', markersize=10, marker='*',
                                      tag='ss21-ss27',
                                      t1='Pre-starburst Phase',
                                      t2='Starburst Phase',
                                      t3='Post-starburst Phase',
                                      t4='Quiescent Phase',
                                      cbarLabelSize=23,
                                      outdir='./',
                                      legendFontSize=18,
                                      saveFig=True)

    # for paper, 3x2 panel for just ss16 and ss22!!!!
    from plot_modules.plot_cloud_prop import plot_stuff_3by2
    ss1, to_plot1, leafdir_out1 = load_pickleTOplot(16)
    ss2, to_plot2, leafdir_out2 = load_pickleTOplot(22)

    fig, ax = plot_stuff_3by2(to_plot1, to_plot2, ls='',
                              markersize=10,
                              marker='*',
                              tag='ss16-ss22',
#                              cmap=cm,
                              sfrlabel=False,
                              cbarLabelSize=23,
                              outdir='./',
                              legendFontSize=18,
                              saveFig=True)
