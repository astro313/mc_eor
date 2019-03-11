'''

plot physical properties of clouds across all snapshots, identified using same way (ncut, steps, etc).

last mod: 18 Aug 2018


'''

import os
from plot_modules.plot_cloud_prop import setup_plot
setup_plot()
import matplotlib.pyplot as plt


def load_pickleTOplot(minss, maxss, fname):
    from plot_modules.plot_cloud_prop import unpack_xy
    from io_modules.leaf_pickle import load_in_pickle_SS

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

    plot_stuff("cloud mass", "mass over jeans mass", leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "jeans mass",  leglabel="sfr: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("cloud mass", "alpha vir", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("size pc", "sigma kms", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("Mach", "SFR young", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("stellar to gas mass", "sigma kms", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
#     plot_stuff("gas sd", "sfr sd", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("gas sd", "sigma kms", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("R2 pc2", "cloud mass", leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff("size pc", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff_3dim("tff Myr", "size pc", "cloud mass", leglabel="ncut: ",
                    to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff('gas sd cgs', 'sigmaSq over size', leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)
    plot_stuff('alpha vir', 'sigmaSq over size', leglabel="ncut: ",
               to_plot=to_plot, outdir=leafdir_out, sfrlabel=True)

    from plot_modules.plot_cloud_prop import get_masses_all_clouds, massFuncUnbinnedCDF, massFuncPDF, massFuncDifferential, CMF
    allmass = get_masses_all_clouds(ss)
    massFuncUnbinnedCDF(allmass, outdir=leafdir_out)
    massFuncPDF(allmass, outdir=leafdir_out)
    massFuncDifferential(allmass, outdir=leafdir_out)
    CMF(allmass, outdir=leafdir_out)


def plot_dis_allSS(minss=16, maxss=28, pattern1='0.32', pattern2='18.96',
             plotdir='./test_brute/leaf_fields_ss16-28/'):

    import numpy as np
    from plot_modules.plot_cloud_prop import get_masses_all_clouds, get_sizes_all_clouds, get_fgas_all_clouds
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    import matplotlib.ticker as ticker

    ss11, to_plot11, _ = load_pickleTOplot(minss, maxss, \
                                           fname=pattern1 + '_10_fields.p')
    ss12, to_plot12, _ = load_pickleTOplot(minss, maxss, \
                                           fname=pattern2 + '_10_fields.p')

    # fig = plt.figure(figsize=(13, 6))
    # ax = fig.subplots(2, 3)

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(13, 3.5))
    fig.subplots_adjust(left=0.12, wspace=0.15, bottom=0.1, top=0.95, hspace=0.0)

    # # separate allmasses into two color scheme
    # med11 = allmass11.median()
    # mask11_low = allmass11 <= med1
    # mask11_high = allmass11 >= med1

    fig.text(0.05, 0.95, r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern1, rotation=90)
    fig.text(0.05, 0.43, r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern2, rotation=90)

    allmass11 = get_masses_all_clouds(ss11)
    ax[0, 0].hist(np.log10(allmass11))
#     ax[0, 0].set_title(r'$n_{\rm cut}$ = %s [cm$^{-3}$]' % pattern1)
    ax[0, 0].set_title(r'$\log M_{\rm cl}$ [M$_{\odot}]$')
    ax[0, 0].set_xlim([5.7, 8.5])
    ax[0, 0].set_xticklabels([])
    allmass12 = get_masses_all_clouds(ss12)
    ax[1, 0].hist(np.log10(allmass12))
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
    allsizes11 = get_sizes_all_clouds(ss11)
    ax[0, 1].hist(allsizes11)
    ax[0, 1].set_xlim([45, 300])
    ax[0, 1].set_xticklabels([])
    ax[0, 1].set_title(r'R [pc]')
    allsizes12 = get_sizes_all_clouds(ss12)
    ax[1, 1].hist(allsizes12)
    ax[1, 1].set_xlim([45, 300])
#    fig.text(0.5, 0.33, r'R [pc]', ha='center', fontsize=20)

    minorLocator = AutoMinorLocator(5)
    ax[1, 0].xaxis.set_minor_locator(minorLocator)
    ax[1, 1].yaxis.set_minor_locator(minorLocator)
    ax[1, 0].xaxis.set_minor_locator(minorLocator)
    ax[1, 1].yaxis.set_minor_locator(minorLocator)

    ax[1, 0].minorticks_on()
    ax[1, 1].minorticks_on()

    # f_gas
    fgas11 = get_fgas_all_clouds(ss11)
    ax[0, 2].hist(fgas11)
    ax[0, 2].set_xlim([0., 1.0])
    ax[0, 2].set_xticklabels([])
    ax[0, 2].set_title(r'$f_{\rm gas}$')
    fgas12 = get_fgas_all_clouds(ss12)
    ax[1, 2].hist(fgas12)
    ax[1, 2].set_xlim([0., 1.0])
#    fig.text(0.51, 0.03, r'$f_{\rm gas}$', ha='center', fontsize=20)

    minorLocator = AutoMinorLocator(5)
    ax[0, 2].xaxis.set_minor_locator(minorLocator)
    ax[1, 2].yaxis.set_minor_locator(minorLocator)
    ax[0, 2].xaxis.set_minor_locator(minorLocator)
    ax[1, 2].yaxis.set_minor_locator(minorLocator)

    ax[0, 2].minorticks_on()
    ax[1, 2].minorticks_on()

    fig.subplots_adjust(hspace=0.45)
    plt.tick_params(which='minor', length=4)

    plt.savefig(plotdir + 'minmaxNcut_basicDistributions.png', \
                bbox_inches='tight')
    plt.close()

    # plot f_gas versus M_gas
    plt.figure()
    plt.plot(allmass11, fgas11, 'o', label=pattern1 + r' cm$^{-3}$')
    plt.plot(allmass12, fgas12, 'D', label=pattern2 + r' cm$^{-3}$')
    plt.legend(loc='best')
    plt.xscale('log')
    plt.xlabel(r"M$_{\rm cl}$ [M$_\odot$]",fontsize=20)
    plt.ylabel(r"f$_{\rm gas}$", fontsize=20)
    plt.tight_layout()
    plt.savefig(plotdir + 'f_gas_M_gas.png', bbox_inches='tight')
    plt.close()


# ---

if __name__ == '__main__':

    # for paper
    plot_dis_allSS()

    # # of MCs w/ mass > 1e8 Msun (for paper)
    # _to, _ = load_pickleTOplot(16, 28, fname="0.32_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="0.53_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="0.88_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="1.47_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="2.45_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="4.08_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="6.81_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="11.36_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count
    # import pdb; pdb.set_trace()

    # _to, _ = load_pickleTOplot(16, 28, fname="18.96_10_fields.p")
    # count = 0
    # for i in _to.iterkeys():
    #     if _to[i]['cloud mass'] > 1.e8:
    #         count += 1
    # print count

    # plotting_procedure(16, 28, pattern="0.32")
    # plotting_procedure(16, 28, pattern="0.53")
    # plotting_procedure(16, 28, pattern="0.88")
    # plotting_procedure(16, 28, pattern="1.47")
    # plotting_procedure(16, 28, pattern="2.45")
    # plotting_procedure(16, 28, pattern="4.08")
    # plotting_procedure(16, 28, pattern="6.81")
    # plotting_procedure(16, 28, pattern="11.36")
    # plotting_procedure(16, 28, pattern="18.96")
