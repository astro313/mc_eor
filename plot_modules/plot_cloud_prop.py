import numpy as np
import matplotlib.pyplot as plt


def setup_plot():
    import matplotlib
    # print(matplotlib.matplotlib_fname())

    matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                                , 'font.size': 10      # points
                                , 'legend.fontsize': 8      # points
                                , 'lines.linewidth': 2       # points
                                , 'axes.linewidth': 1       # points
                                , 'text.usetex': True    # Use LaTeX to layout text
                                , 'font.family': "serif"  # Use serifed fonts
                                , 'xtick.major.size': 10     # length, points
                                , 'xtick.major.width': 1     # points
                                , 'xtick.minor.size': 6     # length, points
                                , 'xtick.minor.width': 1     # points
                                , 'ytick.major.size': 10     # length, points
                                , 'ytick.major.width': 1     # points
                                , 'ytick.minor.size': 6     # length, points

                                , 'ytick.minor.width': 1     # points
                                , 'font.serif': ("Times", "Palatino", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                                })
    return None


def col_f(ii, cm=None):
    if cm is None:
        cm = plt.get_cmap("Magma")
        # cm = plt.get_cmap("cubehelix")
#        cm = plt.get_cmap('gist_heat')
    # cm = plt.get_cmap('gist_rainbow')

    return cm(ii)


def unpack_xy(ss):
    """

    Unpack x, y-axes variables that we want to plot.
    Single snapshot, different ncut.

    Parameters
    ----------
    ss: dict

    Returns
    -------
    to_plot: dict

    """

    pc2kpc = 1.e-3
    cm2km = 1.e-5

    to_plot = {}

    for ks in sorted(ss.iterkeys()):
        _Mj = []
        _MMj = []
        _m = []
        _mach = []
        _mstar = []
        _mstar2mgas = []
        _SFR_young = []
        _SFR_old = []
        _alpha = []
        _mSD = []
        _sizepc = []
        _R2pc2 = []
        _SFRSD = []
        _sigmakms = []
        _tffMyr = []
        _sigmakms = []
        _tffMyr = []
        _mSD_cgs = []
        _sigma2oR = []

        for kkk in ss[ks].iterkeys():
            # print ss[ks][kkk]
            _Mj.append(ss[ks][kkk].M_jeans)
            MMj = ss[ks][kkk].mass_Msun / ss[ks][kkk].M_jeans
            _MMj.append(MMj)
            _mach.append(ss[ks][kkk].Mach)
            _m.append(ss[ks][kkk].mass_Msun)
            _mstar.append(ss[ks][kkk].mstar_Msun_tot)
            _mstar2mgas.append(ss[ks][kkk].s2gR)
            _SFR_young.append(ss[ks][kkk].young_SFR_MsunPyr)
            _SFR_old.append(ss[ks][kkk].old_SFR_MsunPyr)
            _alpha.append(ss[ks][kkk].alpha)
            _mSD.append(ss[ks][kkk].massSD)
            _sizepc.append(ss[ks][kkk].R_pc * 2.0)
            _R2pc2.append((ss[ks][kkk].R_pc * 2.0) ** 2.0)
            _SFRSD.append((ss[ks][kkk].SFR) / 1.0e6 /
                          (np.pi * ss[ks][kkk].R_pc * pc2kpc)**2)
            _sigmakms.append(np.sqrt(ss[ks][kkk].sigmaSq) * cm2km)
            _tffMyr.append(ss[ks][kkk].tff_Myr)
            _mSD_cgs.append(ss[ks][kkk].massSD *
                            ss[ks][kkk].Msun2g / (ss[ks][kkk].pc2cm)**2)
            _sigma2oR.append(ss[ks][kkk].sigmaSq / (1.e5)**2 /
                             (ss[ks][kkk].R_pc * 2.0))

        to_plot[ks] = {}
        to_plot[ks]['cloud mass'] = _m
        to_plot[ks]['Mach'] = _mach
        to_plot[ks]['cloud stellar mass'] = _mstar
        to_plot[ks]['stellar to gas mass'] = _mstar2mgas
        to_plot[ks]['SFR young'] = _SFR_young
        to_plot[ks]['SFR old'] = _SFR_old
        to_plot[ks]['jeans mass'] = _Mj
        to_plot[ks]['mass over jeans mass'] = _MMj
        to_plot[ks]['alpha vir'] = _alpha
        to_plot[ks]['gas sd'] = _mSD
        to_plot[ks]['sfr sd'] = _SFRSD
        to_plot[ks]['size pc'] = _sizepc
        to_plot[ks]['R2 pc2'] = _R2pc2
        to_plot[ks]['sigma kms'] = _sigmakms
        to_plot[ks]['tff Myr'] = _tffMyr
        to_plot[ks]['gas sd cgs'] = _mSD_cgs
        to_plot[ks]['sigmaSq over size'] = _sigma2oR

    return to_plot


def plot_stuff(xstr, ystr, ls='', markersize=7, marker='*',
               leglabel='', tag='', cm='gist_rainbow',
               to_plot=None, save=True, outdir='./',
               sfrlabel=False):
    """

    Parameters
    ----------
    xstr: str
        something like 'cloud mass'

    ystr: str
        something like 'cloud mass', 'alpha vir'

    leglabel: str
         something like "ncut: ", or "snapshot: " for legend

    to_plot: dict

    sfrlabel: bool
        if true, will label marker using SFR instead of snapshot number

    Returns
    -------
    fig, ax

    """
    import pandas as pd
    import glob

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    NUM_COLORS = len(to_plot)
    cm = plt.get_cmap("Magma")

    if xstr == "gas sd" and ystr == "sfr sd":
        _, Heinerman_SigmaGas, _, Heinerman_SigmaSFR = load_Heiderman10()
        ax.plot(Heinerman_SigmaGas, Heinerman_SigmaSFR, marker='.', markersize=7, linestyle='',
                label="MW Heiderman+2010", color='k', alpha=0.8)

        x, y = [10**0.50, 10**4.0], [10**(-2.85), 10**2.1]
        ax.plot(x, y, linestyle='-', color='b',
                linewidth=2, label="Kennicutt 1998")

        # more from high-z literature
        litpath = '/mnt/home/daisyleung/mc_eor/literature/'
        x0901, y0901 = np.loadtxt(
            litpath + "J0901_KS10_points.txt", unpack=True)  # in log
        x14011, y14011 = np.loadtxt(
            litpath + "J14011_KSpoints2.txt", unpack=True)  # not in log
        xrawle, yrawle = np.loadtxt(
            litpath + "Rawle_KSpoints.txt", unpack=True, usecols=(0, 1))  # in log
        xgn20, xgn20err, ygn20, ygn20err = np.loadtxt(
            litpath + "Hodge_resolvedKS.txt", unpack=True)   # not in log
        xegs, xegserr, yegs, yegserr = np.loadtxt(
            litpath + "Genzel_KSpoints.txt", unpack=True)  # in log

        ax.scatter(10**x0901, 10**y0901, label="J0901 @ z=2.26",
                   color='red', marker='o', s=5, facecolors='none', alpha=0.6)
        ax.scatter(x14011, y14011, label="SMM J14011 @ z=2.56",
                   color='darkblue', marker='^', s=13, facecolors='none', alpha=0.8)
        ax.scatter(10**xrawle, 10**yrawle, label="HLS0918 @ z=5.24",
                   color='purple', marker='v', s=10, facecolors='none', alpha=0.8)
        ax.errorbar(xgn20, ygn20, yerr=ygn20err, xerr=xgn20err,
                    label="GN20 @ z=4.05",
                    color='orange', fmt='s', markersize=4.5,
                    markeredgewidth=0.6, mfc='none', elinewidth=0.5)
        ax.errorbar(10**xegs, 10**yegs, yerr=yegserr, xerr=xegserr,
                    label="EGS13011166 @ z=1.53",
                    color='green', fmt='D', markersize=3.5,
                    markeredgewidth=0.6, mfc='none', zorder=0.5, alpha=0.56, elinewidth=0.5)

    if ystr == "alpha vir" and xstr == "cloud mass":
        # Kauffmann+17 Figure 4
        litpath = '/mnt/home/daisyleung/mc_eor/literature/'

        # read data from old Kauffmann paper
        CloudData = pd.read_table(litpath + './Kauffmann17/filter_alpha_vp-paper.dat',
                                  header=0,
                                  delim_whitespace=True)

        # CMZ data: GCMS dendrograms
        DataTable = pd.DataFrame()
        LinewidthSizeFiles = glob.glob(
            litpath + './Kauffmann17/Clumps_*.pickle')
        for File in LinewidthSizeFiles:
            DataTable = DataTable.append(pd.read_pickle(File))

        DataTable['mass'] = 224.46 * (2.0E23 / 1.0E22) * \
            np.pi * (DataTable['r_eff'])**2.
        DataTable['alpha'] = 1.2 * (DataTable['v_rms'])**2. * \
            DataTable['r_eff'] * \
            (DataTable['mass'] / 1.0E3)**(-1)
        DataTable['Mach'] = DataTable['v_rms'] / \
            (0.288 / 2.33**0.5 * (50. / 10.)**0.5)

        np.unique(CloudData['Sample'])

        # show instability range
        ax.fill_between([1.0E-9, 1.0E9],
                        [1.0E-9, 1.0E-9],
                        [2., 2.],
                        facecolor='0.9', edgecolor='none',
                        zorder=0)
        plt.annotate(s='unstable',
                     xy=[2.e6, 0.2],
                     color='0.3', size=22.,
                     ha='left')

        COFilterArray = (CloudData['Sample'] == 'DuvalGRS')
        CloudData.loc[CloudData['Sample'] == 'HeyerGRS'] = 'NaN'

        # plot reference data
        plt.plot(CloudData['Mass'][COFilterArray],
                 CloudData['Alpha'][COFilterArray],
                 'p', markeredgecolor='palegreen', markeredgewidth=2.,
                 markersize=7., markerfacecolor='none', label='Pipe Nebula dense gas')  # Pipe Nebula, C18O (associated with high-mass SF), NH3 (high-density tracer)

        plt.plot(CloudData['Mass'][np.invert(COFilterArray)],
                 CloudData['Alpha'][np.invert(COFilterArray)],
                 'o', color='limegreen', markeredgewidth=0., label='MW clouds dense gas tracers')

        # CMZ data: GCMS dendrograms
        plt.plot(DataTable[DataTable['Target'] != 'SgrD']['mass'],
                 DataTable[DataTable['Target'] != 'SgrD']['alpha'],
                 '+',
                 markersize=13.,
                 markeredgecolor='blue', markeredgewidth=3.,
                 markerfacecolor=(1, 0.7, 0.7),
                 zorder=10, label='CMZ')

        plt.plot(DataTable[DataTable['Target'] == 'SgrD']['mass'],
                 DataTable[DataTable['Target'] == 'SgrD']['alpha'],
                 '+',
                 markersize=13.,
                 markeredgecolor='blue', markeredgewidth=6.,
                 markerfacecolor=(1, 0.7, 0.7),
                 zorder=10, label='')

        plt.plot(DataTable[DataTable['Target'] == 'SgrD']['mass'],
                 DataTable[DataTable['Target'] == 'SgrD']['alpha'],
                 '+',
                 markersize=13.,
                 markeredgecolor='white', markeredgewidth=2.,
                 markerfacecolor=(1, 0.7, 0.7),
                 zorder=10, label='Sgr D outside CMZ')

        ax.annotate(s='"clumps"',
                    xy=[5.0E3, 0.7],
                    color='blue',
                    ha='left')

        # # estimate for dense cores
        # ax.fill_between([1.0E2, 1.0E3],
        #                 1.2 * 0.6**2. * 0.1 /
        #                 (np.array([1.0E2, 1.0E3]) / 1.0E3),
        #                 1.2 * 2.2**2. * 0.1 /
        #                 (np.array([1.0E2, 1.0E3]) / 1.0E3),
        #                 facecolor='blue', edgecolor='none', alpha=0.3,
        #                 zorder=3)
        # plt.plot([1.0E2, 1.0E3],
        #          1.2 * 0.6**2. * 0.1 / (np.array([1.0E2, 1.0E3]) / 1.0E3),
        #          color='blue',
        #          linewidth=3.)
        # annotate(s=r'$\sigma_{\mathdefault{v}} = \mathdefault{0.6 \, km \, s^{-1}}$',
        #          xy=[300., 1. / 1.5 * 1.2 * 0.6**2. * 0.1 / (300. / 1.0E3)],
        #          color='blue',
        #          ha='right')
        # plt.plot([1.0E2, 1.0E3],
        #          1.2 * 2.2**2. * 0.1 / (np.array([1.0E2, 1.0E3]) / 1.0E3),
        #          color='blue',
        #          linewidth=3.)
        # annotate(s=r'$\sigma_{\mathdefault{v}} = \mathdefault{2.2 \, km \, s^{-1}}$',
        #          xy=[100, 1.2 * 1.2 * 2.2**2. * 0.1 / (100. / 1.0E3)],
        #          color='blue',
        #          ha='center')

        # CMZ data: entire clouds
        EntireCloudData = pd.DataFrame()
        EntireCloudData['Target'] = [
            'SgrC', '20kms', '50kms', 'G0.253', 'SgrB1']
        EntireCloudData['Mass'] = [2.5E4, 33.9E4, 6.5E4, 9.3E4, 14.5E4]
        EntireCloudData['Size'] = [1.7, 5.1, 2.7, 2.8, 3.6]
        EntireCloudData['VelocityDispersion'] = [6.5, 10.2, 13.9, 16.4, 13.1]
        EntireCloudData['VirialParameter'] = 1.2 * EntireCloudData['VelocityDispersion']**2. * \
            EntireCloudData['Size'] * \
            (EntireCloudData['Mass'] / 1000.)**-1.
        EntireCloudData['Mach'] = \
            EntireCloudData['VelocityDispersion'] / \
            (0.288 / 2.33**0.5 * (50. / 10.)**0.5)
        EntireCloudData['MeanDensity'] = \
            3.5E4 * (EntireCloudData['Mass'] / 1.0E4) / \
            EntireCloudData['Size']**3.
        EntireCloudData['ThresholdDensitySF'] = \
            EntireCloudData['VirialParameter'] * \
            EntireCloudData['Mach']**2. * \
            EntireCloudData['MeanDensity']

        plt.plot(EntireCloudData['Mass'],
                 EntireCloudData['VirialParameter'],
                 'D',
                 markersize=13.,
                 markeredgecolor='blue', markeredgewidth=3.,
                 markerfacecolor='none',
                 zorder=10, label='CO-based MW')

        plt.annotate(s='entire clouds',
                     xy=[1.0E5, 18.],
                     color='blue',
                     ha='center')

    if xstr == "size pc" and ystr == "sigma kms":
        # Solomon+87: slope=0.5, based on 273 GMCs (superceeded by Heyer+09):
        # y = 0.72 * x**0.5

        # Heyer & Brunt 2004 (27 GMCs in MW): sigma = 0.9 * R^0.56
        x = np.logspace(1, 3, 10)
        yH04 = 0.9 * x**0.56
        ax.plot(x, yH04, linestyle='-.', color='r', linewidth=2,
                label=r'Heyer \& Brunt 2004 $\sigma \propto R^{0.56}$')

        # Bolatto+08: sigma = 0.44 * R^0.6 km/s
        yB08 = 0.44 * x**0.60
        ax.plot(x, yB08, linestyle=':', color='b', linewidth=1.5,
                label=r'Bolatto+08 $\sigma \propto R^{0.60}$')

        # Larson81
        y = 1.10 * x**0.38       # Larson81 Eqn 1
        ax.plot(x, y, color='k', linestyle='--', linewidth=1.5,
                label=r'Larson 1981 $\sigma \propto R^{0.38}$')

        # More data from literature
        litpath = '/mnt/home/daisyleung/mc_eor/literature/'
        xegc, yegc = np.loadtxt(litpath + 'ExtraGalacticGMCs.csv',  # Bolatto08
                                delimiter=',', unpack=True)
        xgc, ygc = np.loadtxt(litpath + 'GalacticCenter.csv',
                              delimiter=',', unpack=True)
        x64, y64 = np.loadtxt(litpath + 'M64.csv', delimiter=',', unpack=True)
        # normalization ~5x higher than found in MW
        xmark, ymark, ymark_err = np.loadtxt(
            litpath + 'SMMJ2135.txt', unpack=True)
        rmark, rmark_err = np.loadtxt(
            litpath + "eyelash_larson.dat", unpack=True, usecols=(5, 6))
        xngc253, yngc253 = np.loadtxt(litpath + 'Leroy15_NGC253.csv',
                                      delimiter=',', unpack=True)

        ax.errorbar(rmark, ymark, yerr=ymark_err, xerr=rmark_err,
                    label="SMM J2135-0102",
                    color='magenta', fmt='D', markersize=3.5,
                    markeredgewidth=1.0)
        ax.scatter(xngc253, yngc253, label="NGC 253",
                   color='red', marker='o', s=7)
        ax.scatter(x64, y64, label="M64", color='orange',
                   marker='^', s=10)
        ax.scatter(xgc, ygc, label="Heyer Galactic Center",
                   color='b', marker='o', s=7)
        ax.scatter(xegc, yegc, label="Bolatto+08: Extra-Galactic GMCs",
                   color='k', marker='.', s=10)

    if xstr == "gas sd cgs" and ystr == "sigmaSq over size":
        # Heyer+09 GRS data
        litpath = '/mnt/home/daisyleung/mc_eor/literature/'
        xxx, yyy = np.loadtxt(litpath + "GRS.txt", unpack=True)
        ax.scatter(xxx, yyy, marker='.', color='k', s=7, label='Heyer+09 GRS')

        POverKb = [1.e4, 1.e5, 1.e6, 1.e7]
        sd_cgs = np.logspace(-3, 3.0, 30)
        Plabel = ['4', '5', '6', '7']
        cm2pc = 1. / 3.086E+18

        def V0_sq_func(pressure, sd_cgs, Gamma=3. / 5):
            """

            Gamma = 3./5 for uniform density sphere

            V_0^2 = sigma^2/R = 1/3 (pi * Gamma * G * Sigma + 4 * P / Sigma)
            K cm^-3, Elmegreen+89: typical Pressure for neutral ISM ~9 x 1E3 K
            cm^-3

            """
            import astropy.constants as C

            print np.pi * Gamma * C.G.cgs.value * sd_cgs
            print pressure / sd_cgs
            import pdb
            pdb.set_trace()

            V0_sq = 1. / 3 * (np.pi * Gamma * C.G.cgs.value * sd_cgs +
                              4. * pressure / sd_cgs)
            print V0_sq
            return V0_sq

        # for ip, Pext in enumerate(POverKb):
        #     V0_sq = V0_sq_func(Pext, sd_cgs)
        #     ax.plot(sd_cgs, V0_sq/(1.e5)**2 / cm2pc,
        #             label=r'Log P = ' + Plabel[ip] + r"K cm$^{-2}$")

    # --- my clouds ----
    ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS)
                                for i in range(NUM_COLORS)])

    if sfrlabel:
        # load in SFR of each snapshot (averaged over 4 Myr, see load_misc.py)
        from io_modules.load_misc import load_SFR_perSS
        sfr = load_SFR_perSS()

    for ks in iter(sorted(to_plot.iterkeys())):
        _x = to_plot[ks][xstr]
        _y = to_plot[ks][ystr]
        if sfrlabel:
            ax.plot(_x, _y, ls=ls, markersize=markersize,
                    marker=marker,
                    label="SFR: " + "{0:d}".format(int(sfr[ks])))
        else:
            ax.plot(_x, _y, ls=ls, markersize=markersize,
                    marker=marker, label=leglabel + ks)

    if(xstr == 'cloud mass'):
        ax.set_xscale("log")
        ax.set_xlabel(r"$M_{\rm cl}$ [M$_{\odot}$]")
        ax.set_xlim(1.0e3, 1.0e8)

    if(xstr == 'cloud stellar mass'):
        ax.set_xscale("log")
        ax.set_xlabel(r"$M_{\rm cl}^*$ [M$_{\odot}$]")
        ax.set_xlim(1.0e5, 1.0e10)

    if xstr == "Mach":
        ax.set_xlabel("Mach ")

    if xstr == "stellar to gas mass":
        # ax.set_xlim(0.0, 0.5)
        ax.set_xlabel(r"$M_*/M_{\rm gas}$")
        ax.set_xscale("log")

    if(xstr == "gas sd"):
        ax.set_xscale("log")
        ax.set_xlim(1.e2, 1.e3)
        ax.set_xlabel(r"$\Sigma_{\rm gas}$ [M$_{\odot}$ pc$^{-2}$]")

    if xstr == 'size pc':
        ax.set_xscale("log")
        ax.set_xlabel("Cloud Size [pc]")
        ax.set_xlim(1.0, 2.e3)

    if xstr == 'R2 pc2':
        ax.set_xscale("log")
        ax.set_xlim(10.0, 1.e8)
        ax.set_xlabel(r"$R^2$ [pc$^{2}$]")

    if xstr == 'gas sd cgs':
        ax.set_xscale("log")
        ax.set_xlim(10**-3, 10**1.0)
        ax.set_xlabel(r"$\Sigma_{\rm gas}$ [g cm$^{-2}$]")

    if ystr == 'sigma kms':
        ax.set_yscale("log")
        ax.set_ylabel(r"$\sigma$ [km s$^{-1}$]")
        ax.set_ylim(0.1, 7.e2)

    if ystr == "SFR young":
        ax.set_yscale("log")
        ax.set_ylabel("SFR from existing young stars [Msun " + r"yr$^{-1}$]")

    if(ystr == "gas sd"):
        ax.set_ylim(1.0, 1.e4)

    if ystr == "sfr sd":
        ax.set_yscale("log")
        ax.set_ylabel(
            r"$\Sigma_{\rm SFR}$ [M$_{\odot}$ yr$^{-1}$ kpc$^2$]")

    if(ystr == 'cloud mass'):
        ax.set_yscale("log")
        ax.set_ylabel(r"$M_{\rm cl}$ [M$_{\odot}$]")
        ax.set_ylim(1.0e3, 1.0e8)

    if(ystr == 'mass over jeans mass'):
        ax.set_yscale("log")
        ax.set_ylabel(r"$M_{\rm cl} / $M$_J$")
#        ax.set_ylim(10.0, 5.e6)

    if(ystr == 'alpha vir'):
        ax.set_yscale("log")
        ax.set_ylabel(r"$\alpha_{\rm vir}$")
        ax.set_ylim(0.02, 2.e2)

        if(xstr == 'cloud mass'):
            ax.set_xlim(0.02, 1.0e8)    # to accomodate Kauffmann+17 data

    if ystr == 'jeans mass':
        ax.set_yscale('log')

    if ystr == 'sigmaSq over size':
        ax.set_yscale("log")
        ax.set_ylim(10**-1.5, 10**2.5)
        ax.set_ylabel(r"$\sigma^2/R$ [km$^2$ s$^{-2}$ pc$^{-1}$]")

    if leglabel is not '':

        def sort_f(t, ss, leglabel='ncut'):
            # ss = {'CO-based MW': 1, 'Sgr D outside CMZ': 2,
            #       'MW clouds dense gas tracers': 3, 'Pipe Nebula dense gas': 4}
            try:
                out = ss[t]
            except KeyError:
                print t
                val = t[t.index(leglabel) + len(leglabel)+1:]
                print val
                val = float(val)
                out = int(100 * (1.+ val))
            return out

        handles, labels = ax.get_legend_handles_labels()
        order = range(len(labels))

        if xstr == "cloud mass" and ystr == "alpha vir":
            ss = {'CMZ': 1,
                  'CO-based MW': 2,
                  'Sgr D outside CMZ': 3,
                  'MW clouds dense gas tracers': 4,
                  'Pipe Nebula dense gas': 5}
            order = sorted(order, key=lambda t: sort_f(labels[t], ss, leglabel))

        if xstr == "size pc" and ystr == "sigma kms":
                # Shrink current axis by 20%
            box = ax.get_position()
            # ax.set_position([box.x0, box.y0 + 0.25, box.width, box.height])
            # ax.legend(loc="upper center", ncol=4, fontsize=9,
            #           bbox_to_anchor=(0.5, -0.18))
            ax.legend([handles[idx] for idx in order],
                      [labels[idx] for idx in order],
                      loc="upper center", ncol=4,
                      fontsize=9,
                      bbox_to_anchor=(0.5, -0.18))
        elif xstr == "cloud mass" and ystr == 'alpha vir':
            ax.legend([handles[idx] for idx in order],
                      [labels[idx] for idx in order], loc='best', ncol=2, fontsize=10)
        else:
            ax.legend([handles[idx] for idx in order],
                      [labels[idx] for idx in order], loc='best', fontsize=10)

    ax.tick_params(axis='both', which='both')   # direction='in'
    plt.tight_layout()

    if save:
        name_out = ystr.replace(' ', '-') + '_' + \
            xstr.replace(' ', '-')
        fig.savefig(outdir + name_out + tag + '.png', bbox_inches="tight")
    else:
        plt.show()

    return fig, ax


def plot_stuff_3dim(xstr, ystr, zstr, ls='', markersize=7, marker='*',
                    leglabel='', tag='', cm='gist_rainbow',
                    to_plot=None, save=True, outdir='./', sfrlabel=False):
    from matplotlib.ticker import NullFormatter

    if sfrlabel:
        # load in SFR of each snapshot (averaged over 4 Myr, see load_misc.py)
        from io_modules.load_misc import load_SFR_perSS
        sfr = load_SFR_perSS()

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    NUM_COLORS = len(to_plot)
    cm = plt.get_cmap(cm)

    # --- my clouds ----
    ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS)
                                for i in range(NUM_COLORS)])

    for ks in iter(sorted(to_plot.iterkeys())):
        _x = to_plot[ks][xstr]
        _y = to_plot[ks][ystr]
        _z = to_plot[ks][zstr]

        if sfrlabel:
            cax = ax.scatter(_x, _y, s=np.array(_z) / np.array(_z).min(),
                             marker=marker,
                             label="SFR: " + "{0:d}".format(int(sfr[ks])))
        else:
            cax = ax.scatter(_x, _y, s=np.array(_z) / np.array(_z).min(),
                             marker=marker, label=leglabel + ks)

    if(xstr == 'tff Myr'):
        ax.set_xlabel(r"$t_{\rm ff}$ [Myr]")
        ax.set_xlim(1.0, 5.0)

    if ystr == "size pc":
        ax.set_yscale("log")
        ax.set_ylim(10.0, 1.e4)
        ax.set_ylabel(r"$R$ [pc]")
        # ax.yaxis.set_major_locator(plt.MaxNLocator(1))
        # ax.yaxis.set_minor_locator(plt.MaxNLocator(1))
        ax.yaxis.set_minor_formatter(NullFormatter())

    if leglabel is not '':
        ax.legend(loc='best')

    ax.tick_params(axis='both', which='both')   # direction='in'
    plt.tight_layout()

    if save:
        name_out = ystr.replace(' ', '-') + '_' + \
            xstr.replace(' ', '-') + '_' + \
            zstr.replace(' ', '-') + '_'
        fig.savefig(outdir + name_out + tag + '.png', bbox_inches="tight")
    else:
        plt.show()

    return fig, ax, cax


def load_Heiderman10():
    """ Read in data from table 1 of Heiderman+10 regarding the MW. """

    # Table 1 in Heiderman 2010
    cloud_name = ["Cha II", "lup I", 'lup II', 'lup IV',
                  'Oph', 'Per', 'Ser',
                  'AurN', 'Aur', 'Cep',
                  'Cha III', 'Cha I',
                  'CrA', 'IC5146E', 'IC5146NW',
                  'Lup VI', 'LupV', 'Mus',
                  'Sco', 'Ser-Aqu']

    Heinerman_units = ["Sigma Gas = [Msun / pc^2]",
                       "SFR = [Msun / Myr]",
                       "Sigma SFR = [Msun / yr kpc^2]"]

    Heinerman_mass = [637, 513, 912, 189, 3120,
                      6590, 2340, 224, 4620, 2610,
                      1330, 857, 279, 3370,
                      5180, 455, 705, 335, 621, 24400]

    Heinerman_SigmaGas = [64.3, 57.9, 59.2, 75.0,
                          105, 90.0, 138, 92.9, 92.4,
                          68.7, 47.5, 91.1, 92.1, 54.9,
                          59.1, 67.5, 60.3, 49.1, 85.2, 136]

    Heinerman_SFR = [6.00, 3.25, 17.0, 3.0, 72.5, 96.2,
                     56, 0.5, 42.7, 29.5, 1.0, 22.2, 10.2, 23.2,
                     9.50, 11.2, 10.7, 3.0, 2.5, 360]

    Heinerman_SigmaSFR = [0.605, 0.367, 1.10, 1.19, 2.45,
                          1.31, 3.29, 0.207, 0.854, 0.776, 0.0357,
                          2.36, 3.37, 0.378, 0.108, 1.66, 0.915,
                          0.440, 0.343, 2.01]

    return Heinerman_mass, Heinerman_SigmaGas, Heinerman_SFR, Heinerman_SigmaSFR


def mass_function(data, logged=False,
                  nbins=35, verbose=True):
    """
    Calculates a mass function from data.

    Parameters
    ----------
    data: array

    logged: boolean
        if False, will take log inside function

    nbins: int

    verbose: boolean


    Return
    ------
    mf / dm:
        dN/dlnM, differential mass function
    mbin:
        mass bins

    """

    nbins = int(nbins)

    # if log have been taken from the data or not
    if not logged:
        dat = np.log(data)
    else:
        dat = data
    del data

    # number of object
    nobj = len(dat)

    # bins
    mmax = dat.max()
    mmin = dat.min()
    dm = (mmax - mmin) / float(nbins)
    # mid points
    mbin = (np.arange(nbins) + 0.5) * dm + mmin

    # count up mass in each bin
    mf = np.zeros(nbins)
    # which bin each data point belongs to
    ibin = np.floor((dat - mmin) / dm)
    # make a mask to keep bins given user-defind nbins
    mask = (ibin >= 0) & (ibin < nbins)
    # sum up each bin
    wght = np.ones(len(ibin[mask]))
    for i in range(nbins):
        mf[i] = np.sum(wght[np.array(ibin[mask], dtype=int) == i])
    mf = np.array(mf)

    if not logged:
        mbin = np.e ** mbin

    if verbose:
        print 'Number of object = %i' % nobj
        print 'dlnM =', dm
        print 'min = %f, max = %f' % (mmin, mmax)
        print 'Results:\n', mbin, mf / dm

    return mbin, mf / dm


def get_masses_all_clouds(ss):
    """
    For plotting cloud mass distribution.
    Gather all clouds in ss.

    Parameters
    ----------
    ss: dict

    Returns
    -------
    allmasses: array
        cloud masses

    """

    allmasses = []
    for snap in ss.iterkeys():
        for snapleafs in ss[snap].iterkeys():
            allmasses.append(ss[snap][snapleafs].mass_Msun)
    return np.array(allmasses)


def massFuncUnbinnedCDF(allmasses, save=True, outdir='./', tag=''):
    """ unbinned CDF """

    X2 = np.sort(allmasses)
    F2 = np.array(range(len(allmasses))) / float(len(allmasses))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X2, F2, label='unbinned CDF', lw=2, alpha=1, zorder=2)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$M_{\rm cl}$ [M$_\odot$]")
    ax.set_ylabel("CDF")

    ax.set_xlim(allmasses.min(), allmasses.max())

    ax.legend(loc="best")
    plt.tight_layout()
    if save:
        fig.savefig(outdir + 'MassDistribution_' +
                    tag + '.png', bbox_inches="tight")
    else:
        plt.show()

    return None


def massFuncPDF(allmasses, nbins=200, normed=True, save=True, outdir='./', tag=''):
    """ PDF of all cloud masses """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(allmasses, bins=nbins, normed=normed)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(r"$M_{\rm cl}$ [M$_\odot$]")
    ax.set_ylabel("PDF")

    ax.set_xlim(allmasses.min(), allmasses.max())

    ax.legend(loc="best")
    plt.tight_layout()
    if save:
        fig.savefig(outdir + 'MassDistributionPDF_' + tag + '.png',
                    bbox_inches="tight")
    else:
        plt.show()


def massFuncDifferential(allmasses, logged=False, nbins=8, save=True, outdir='./', tag='', verbose=False):
    """

    Plot dN/dln M.

    Parameters
    ----------
    allmasses: array

    logged:
        whether allmasses is already in log-scale

    nbins: int


    """

    fig = plt.figure()
    ax = fig.add_subplot(121)
    # allmasses_low = np.array(allmasses[allmasses < 5.e6])
    # allmasses_low_log = np.log(allmasses_low)

    # h = 1.06 * np.std(allmasses_low) * len(allmasses_low)**-0.2
    # nbin = int((allmasses_low.max() - allmasses_low.min())/h)
    # print nbin
    # nbin = 35

    # mbin, df = mass_function(allmasses_low_log, logged=True, nbins=nbin)
    # ax.plot(mbin, df)

    # since it's pretty bimodal?!
    mbin, df = mass_function(allmasses[allmasses < 2.e7],
                             logged=logged,
                             nbins=nbins)
    ax.plot(mbin / 1.e7, df)
    ax.set_xlabel(r"$M_{\rm cl}$ [$\times$ 10$^7$ M$_\odot$]")
    ax.set_ylabel("dN/d ln M")

    ax = fig.add_subplot(122)
    mbin, df = mass_function(allmasses[allmasses > 2.e7],
                             logged=logged,
                             nbins=8,
                             verbose=verbose)
    ax.plot(mbin / 1.e7, df)

    # ax.set_xscale("log")
    # ax.set_yscale("log")

    # ax.set_xlim(allmasses.min(), allmasses.max())
    plt.tight_layout()
    if save:
        fig.savefig(outdir + 'MassDistribution_differential_' + tag + '.png',
                    bbox_inches="tight")
    else:
        plt.show()
