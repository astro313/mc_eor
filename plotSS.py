'''

plot physical properties of clouds across all snapshots, identified using same way (ncut, steps, etc).

last mod: 18 Aug 2018


'''

def load_pickleTOplot(minss, maxss, fname):
    from plot_modules.plot_cloud_prop import unpack_xy
    from io_modules.leaf_pickle import load_in_pickle_SS
    import os

    here = os.path.dirname(os.path.abspath(__file__))
    pickle_base_fold = here + '/test_brute/'

    ss = load_in_pickle_SS(pickle_base_fold + 'leaf_fields_', \
                                fname, minss, maxss)
    to_plot = unpack_xy(ss)

    outleafdir = pickle_base_fold + 'leaf_fields_ss' + \
        str(minss) + '-' + str(maxss) + '/'
    leafdir_out = outleafdir + fname[:fname.find('.p')] + '/'
    if not os.path.isdir(leafdir_out):
        os.system('mkdir -p ' + leafdir_out)

    return ss, to_plot, leafdir_out


def plotting_procedure(minss, maxss, pattern):
    from plot_modules.plot_cloud_prop import plot_stuff, plot_stuff_3dim

    fname = pattern + '_10_fields.p'
    ss, to_plot, leafdir_out = load_pickleTOplot(minss, maxss, fname)

    plot_stuff("stellar to gas mass", "cloud mass", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("stellar to gas mass", "alpha vir", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)

    plot_stuff("alpha vir", "SFR young", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "mass over jeans mass", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "jeans mass",  leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "alpha vir", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("size pc", "sigma kms", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("Mach", "SFR young", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("Mach", "SFR old", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("stellar to gas mass", "sigma kms", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("gas sd", "sfr sd", leglabel="sfr: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("gas sd", "sigma kms", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("R2 pc2", "cloud mass", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("size pc", "cloud mass", leglabel="sfr: ",
                to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "Mach", leglabel="sfr: ",
                to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("size pc", "weighted density", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out,  sfrlabel=True)
    plot_stuff('sigma kms NT', 'sigma kms bulk', leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)

    plot_stuff_3dim("tff Myr", "size pc", "cloud mass", leglabel="sfr: ",
                    to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff('gas sd cgs', 'sigmaSq over size', leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff('alpha vir', 'sigmaSq over size', leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)


    # from plot_modules.plot_cloud_prop import get_masses_all_clouds, massFuncUnbinnedCDF, massFuncPDF, massFuncDifferential, CMF
    # allmass = get_masses_all_clouds(ss)
    # massFuncUnbinnedCDF(allmass, outdir=leafdir_out)
    # massFuncPDF(allmass, outdir=leafdir_out)
    # massFuncDifferential(allmass, outdir=leafdir_out)
    # CMF(allmass, outdir=leafdir_out)


def plot_dis_allSS(minss=16, maxss=28, pattern1='0.32', pattern2='18.96',
             plotdir='./test_brute/leaf_fields_ss16-28/'):

    import numpy as np
    from plot_modules.plot_cloud_prop import get_masses_all_clouds, get_sizes_all_clouds, get_fgas_all_clouds
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    import matplotlib.ticker as ticker
    import matplotlib.pyplot as plt

    ss11, to_plot11, _ = load_pickleTOplot(minss, maxss, \
                                           fname=pattern1 + '_10_fields.p')
    ss12, to_plot12, _ = load_pickleTOplot(minss, maxss, \
                                           fname=pattern2 + '_10_fields.p')

    # fig = plt.figure(figsize=(13, 6))
    # ax = fig.subplots(2, 3)

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(16, 5))
    fig.subplots_adjust(left=0.12, wspace=0.18, bottom=0.1, top=0.95, hspace=0.2)

    # # separate allmasses into two color scheme
    # med11 = allmass11.median()
    # mask11_low = allmass11 <= med1
    # mask11_high = allmass11 >= med1

    fig.text(0.05, 0.95, r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern1, rotation=90, fontsize=20)
    fig.text(0.05, 0.43, r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern2, rotation=90, fontsize=20)

    binwidth = 0.1
    allmass11 = np.log10(get_masses_all_clouds(ss11))
    weights = np.ones_like(allmass11)/float(len(allmass11))
    ax[0, 0].hist(allmass11, bins=np.arange(min(allmass11), max(allmass11) +binwidth, binwidth), weights=weights)
#     ax[0, 0].set_title(r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern1)
    ax[0, 0].set_title(r'$\log M_{\rm cl}$ [M$_{\odot}]$', fontsize=20)
    ax[0, 0].set_ylabel(r'$N$', fontsize=20)
    ax[0, 0].set_xlim([5.7, 8.5])
    ax[0, 0].set_xticklabels([])
    allmass12 = np.log10(get_masses_all_clouds(ss12))
    weights = np.ones_like(allmass12)/float(len(allmass12))
    ax[1, 0].hist(allmass12, bins=np.arange(min(allmass12), max(allmass12) +binwidth, binwidth), weights=weights)
    ax[1, 0].set_ylabel(r'$N$', fontsize=20)
    ax[1, 0].set_xlim([5.7, 8.5])
#     ax[0, 1].set_title(r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern2)
#    fig.text(0.51, 0.62, r'$\log $ M$_{\rm cl}$ [M$_{\odot}]$', ha='center', fontsize=20)

    minorLocator = AutoMinorLocator(5)
    ax[0, 0].xaxis.set_minor_locator(minorLocator)
    ax[0, 0].yaxis.set_minor_locator(minorLocator)
    ax[0, 1].xaxis.set_minor_locator(minorLocator)
    ax[0, 1].yaxis.set_minor_locator(minorLocator)

    ax[0, 0].minorticks_on()
    ax[0, 1].minorticks_on()

    # size
    binwidth = 10
    allsizes11 = get_sizes_all_clouds(ss11)
    weights = np.ones_like(allsizes11)/float(len(allsizes11))
    ax[0, 1].hist(allsizes11, bins=np.arange(min(allsizes11), max(allsizes11) +binwidth, binwidth), weights=weights)
    ax[0, 1].set_xlim([0, 260])
    ax[0, 1].set_xticklabels([])
    ax[0, 1].set_title(r'R [pc]', fontsize=20)
    allsizes12 = get_sizes_all_clouds(ss12)
    weights = np.ones_like(allsizes12)/float(len(allsizes12))
    ax[1, 1].hist(allsizes12, bins=np.arange(min(allsizes12), max(allsizes12) +binwidth, binwidth), weights=weights)
    ax[1, 1].set_xlim([0, 260])
#    fig.text(0.5, 0.33, r'R [pc]', ha='center', fontsize=20)

    minorLocator = AutoMinorLocator(5)
    ax[1, 0].xaxis.set_minor_locator(minorLocator)
    ax[1, 0].yaxis.set_minor_locator(minorLocator)
    ax[1, 1].xaxis.set_minor_locator(minorLocator)
    ax[1, 1].yaxis.set_minor_locator(minorLocator)

    ax[1, 0].minorticks_on()
    ax[1, 1].minorticks_on()

    # to only use integer yaxis tick
    from matplotlib.ticker import MaxNLocator
    ax[1, 1].yaxis.set_major_locator(MaxNLocator(2, integer=True))

    # f_gas
    binwidth = 0.1
    fgas11 = get_fgas_all_clouds(ss11)
    weights = np.ones_like(fgas11)/float(len(fgas11))
    ax[0, 2].hist(fgas11, bins=np.arange(min(fgas11), max(fgas11) +binwidth, binwidth),  weights=weights)
    ax[0, 2].set_xlim([0., 1.0])
    ax[0, 2].set_xticklabels([])
    ax[0, 2].set_title(r'$f_{\rm gas}$', fontsize=20)
    fgas12 = get_fgas_all_clouds(ss12)
    weights = np.ones_like(fgas12)/float(len(fgas12))
    ax[1, 2].hist(fgas12, bins=np.arange(min(fgas12), max(fgas12) +binwidth, binwidth),  weights=weights)
    ax[1, 2].set_xlim([0., 1.0])
#    fig.text(0.51, 0.03, r'$f_{\rm gas}$', ha='center', fontsize=20)

    minorLocator = AutoMinorLocator(5)
    ax[0, 2].xaxis.set_minor_locator(minorLocator)
    ax[0, 2].yaxis.set_minor_locator(minorLocator)
    ax[1, 2].xaxis.set_minor_locator(minorLocator)
    ax[1, 2].yaxis.set_minor_locator(minorLocator)

    ax[0, 2].minorticks_on()
    ax[1, 2].minorticks_on()

    fig.subplots_adjust(hspace=0.45)
    plt.tick_params(which='minor', length=4)

    plt.savefig(plotdir + 'minmaxNcut_basicDistributions.pdf', \
                bbox_inches='tight')
    plt.close()

    # plot f_gas versus M_gas
    plt.figure()
    plt.plot(allmass11, fgas11, 'o', label=pattern1 + r' cm$^{-3}$')
    plt.plot(allmass12, fgas12, 'D', label=pattern2 + r' cm$^{-3}$')
    plt.legend(loc='best')
    plt.xlabel(r"M$_{\rm cl}$ [M$_\odot$]",fontsize=20)
    plt.ylabel(r"f$_{\rm gas}$", fontsize=20)
    plt.tight_layout()
    plt.savefig(plotdir + 'f_gas_M_gas.png', bbox_inches='tight')
    plt.close()


# ---

if __name__ == '__main__':

    import matplotlib
    matplotlib.use('agg')

    import numpy as np

    import matplotlib.pyplot as plt
    from plot_modules.plot_cloud_prop import setup_plot, get_P_NT_all_clouds
    cm = setup_plot()

    # for paper
    # minss=16
    # maxss=28
    # pattern1='0.32'
    # pattern2='18.96'
    # plotdir='./test_brute/leaf_fields_ss16-28/'
    # plot_dis_allSS(minss, maxss, pattern1, pattern2, plotdir)
    # ss11, to_plot11, _ = load_pickleTOplot(minss, maxss, \
    #                                        fname=pattern1 + '_10_fields.p')
    # ss12, to_plot12, _ = load_pickleTOplot(minss, maxss, \
    #                                        fname=pattern2 + '_10_fields.p')
    # print(get_P_NT_all_clouds(ss11).min(), get_P_NT_all_clouds(ss11).max(), np.median(get_P_NT_all_clouds(ss11)))
    # print(get_P_NT_all_clouds(ss12).min(), get_P_NT_all_clouds(ss12).max(), np.median(get_P_NT_all_clouds(ss12)))

    # # of MCs w/ mass > 1e8 Msun (for paper)
    # _, _to, _ = load_pickleTOplot(16, 28, fname="0.32_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="0.53_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="0.88_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="1.47_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count       # 0 MCCs w/ Mcl>1e8 Msun with this ncut already
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="2.45_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="4.08_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="6.81_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="11.36_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _, _to, _ = load_pickleTOplot(16, 28, fname="18.96_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     for ii in range(len(_to[i]['cloud mass'])):
    #         if _to[i]['cloud mass'][ii] > 1.e8:
    #             count += 1
    # print count
    # import pdb; pdb.set_trace()

    plotting_procedure(16, 28, pattern="0.32")
    plotting_procedure(16, 28, pattern="0.53")
    plotting_procedure(16, 28, pattern="0.88")
    plotting_procedure(16, 28, pattern="1.47")
    plotting_procedure(16, 28, pattern="2.45")
    plotting_procedure(16, 28, pattern="4.08")
    plotting_procedure(16, 28, pattern="6.81")
    plotting_procedure(16, 28, pattern="11.36")
    plotting_procedure(16, 28, pattern="18.96")

    # for paper, 3x2 panel for just
    plt.close('all')

    from plot_modules.plot_cloud_prop import setup_plot, plot_stuff, plot_stuff_3by2

    fname1 = '0.53_10_fields.p'
    fname2 = '18.96_10_fields.p'
    _, to_plot1, leafdir_out1 = load_pickleTOplot(16, 28, fname1)
    _, to_plot2, leafdir_out2 = load_pickleTOplot(16, 28, fname2)

    bubu = [fname1[:fname1.find('_')], fname2[:fname2.find('_')]]

    fig, ax = plot_stuff_3by2(to_plot1, to_plot2, ls='',
                              markersize=10,
                              marker='*',
                              tag='allss',
#                              cmap=cm,
                              sfrlabel=bubu,
                              cbarLabelSize=25,
                              outdir='./',
                              legendFontSize=20,
                              axwidth=1.5,
                              tickwidth=2,
                              labelsize=25,
                              ticklabsize=20,
                              saveFig=True)