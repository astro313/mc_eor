'''

plot physical properties of clouds across all snapshots, identified using same way (ncut, steps, etc).

last mod: 15 July 2018

NOTE
----
- hard-code to get dx for now...

'''

import matplotlib
print(matplotlib.matplotlib_fname())

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
import matplotlib.pyplot as plt
import cPickle as pickle
from snapshot28_leafprop import Cloud
import numpy as np
import os
from fetch_gal_fields import get_units
import pymses

# defined as ... '{0:.2f}'.format(args.ncut) + '_' + str(args.step) + '_'
# + str(args.Nmin) from snapshot28_leafprop.py
# fname = "0.03_5_3_fields.p"
# fname = "0.1_5_10_fields.p"
fname = "0.30_5_10_fields.p"
# fname = "1.00_5_10_fields.p"
# fname = "1.00_3_10_fields.p"

outleafdir = 'leaf_fields_ss16-28/'
leafdir_out = outleafdir + fname[:fname.find('.p')] + '/'
if not os.path.isdir(leafdir_out):
    os.system('mkdir -ip ' + leafdir_out)


def get_dx(ssnum):
    # hard-code to get dx for now...
    # saved in fetch_gal_fields.py
    ro = pymses.RamsesOutput("output", ssnum)
    ds = np.load('snapshot'+ str(ssnum) +'_center_fields0123456-15.npz')
    dx_vector = ds['dx_vector']

    originalLevels = np.log2(1. / np.unique(dx_vector))
    highestRes = 2.**(originalLevels.max() * -1)
    # after finely sample unigrid
    dx_pc = highestRes * get_units(ro=ro)['dx'][0]
    return dx_pc

ss = {}
for issnum in range(16, 29):
    leafdir = 'leaf_fields_' + str(issnum) + '/'
    try:
        leaf_fields = pickle.load(open(leafdir + fname, "rb"))
    except:
        continue

    Cloud_dict = {}
    for kkk in leaf_fields.iterkeys():
        Cloud_dict[kkk] = Cloud(get_dx(issnum), leaf_fields[kkk], int(kkk))
        print Cloud_dict[kkk]

    ss[str(issnum)] = Cloud_dict


NUM_COLORS = len(range(16, 29))
cm = plt.get_cmap('gist_rainbow')

# plotting
plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        print ss[ks][kkk]
        MMj = ss[ks][kkk].mass_Msun / ss[ks][kkk].M_jeans
        _y.append(MMj)
        _x.append(ss[ks][kkk].mass_Msun)
    ax.plot(_x, _y, ls='', markersize=7, marker='*', label="snapshot: " + ks)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$M_{\rm cl}$ [M$_{\odot}$]")
ax.set_ylabel(r"$M_{\rm cl} / $M$_J$")

ax.set_xlim(1.0e3, 1.0e8)
ax.set_ylim(10.0, 5.e6)

ax.legend(loc='best')

plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'my_cloud_sample.png', bbox_inches="tight")


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        _y.append(ss[ks][kkk].alpha)
        _x.append(ss[ks][kkk].mass_Msun)
    ax.plot(_x, _y, ls='', markersize=7, marker='*', label="snapshot: " + ks)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Cloud Mass [M$_{\odot}$]")
ax.set_ylabel(r"$\alpha_{\rm vir}$")

ax.set_xlim(1.e3, 1.e8)
ax.set_ylim(0.5, 1.e2)

ax.legend(loc="best")

plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + "alphaVir_CloudMass.png", bbox_inches="tight")

# ---
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    # _y2 = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].mass_Msun)
        _y.append(np.array(ss[ks][kkk].SFR) / 1.0e6)  # to be replaced using LIR
        # _y2.append(np.array(ss[ks][kkk].SFR_JML) / 1.0e6)

    ax.plot(_x, _y, ls='', markersize=7, marker='*', label="snapshot: " + ks)
    # ax.plot(_x, _y2, ls='', markersize=8, marker="D") #, label="SFR JML06")


ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlim(1.e3, 1.e8)
ax.set_ylim(1.e-3, 10.0)

ax.set_xlabel(r"Cloud Mass [M$_{\odot}$]")
ax.set_ylabel(r"SFR [M$_{\odot}$ yr$^{-1}$]")

ax.legend(loc='best')

# ax.set_xlim(1.0e3, 1.0e7)
# ax.set_ylim(1.0e-5, 1.0e1)
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + "SFR_CloudMass.png", bbox_inches="tight")


fig = plt.figure()
ax = fig.add_subplot(111)

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


x, y = [10**0.50, 10**4.0], [10**(-2.85), 10**2.1]

pc2kpc = 1.e-3


ax.plot(x, y, '-b', linewidth=2, label="Kennicutt 1998")  # \citep{SK}
ax.plot(Heinerman_SigmaGas, Heinerman_SigmaSFR, 'b.',
        linewidth=2, label="MW Heiderman+2010", markersize=10)

ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].massSD)
        _y.append((ss[ks][kkk].SFR) / 1.0e6 /
                  (np.pi * ss[ks][kkk].R_pc * pc2kpc)**2)

    ax.plot(_x, _y, ls='',
            markersize=8, marker='*', label='snapshot: ' + str(ks))


ax.set_xscale("log")
ax.set_yscale("log")

ax.tick_params(axis='both', which='both')   # direction='in'

ax.set_xlabel(r"$\Sigma_{\rm gas}$ [M$_{\odot}$ pc$^{-2}$]")
ax.set_ylabel(
    r"$\Sigma_{\rm SFR}$ [M$_{\odot}$ yr$^{-1}$ kpc$^2$]")

ax.legend(loc="best")
# ax.set_xlim(5.0e0, 1.0e3)
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'SurfSFR_SurfGas.png', bbox_inches="tight")


# --- Larson's relations ---
fig = plt.figure()
ax = fig.add_subplot(111)

cm2km = 1.e-5

x = np.logspace(1, 3, 10)
#y = 1.3e-2*x**0.38
y = 1.10 * x**0.38       # Larson81 Eqn 1, where x is in units of pc

# Solomon+87: slope=0.5, based on 273 GMCs (superceeded by Heyer+09):
# y = 0.72 * x**0.5


# Heyer & Brunt 2004 (27 GMCs in MW): sigma = 0.9 * R^0.56
yH04 = 0.9 * x**0.56
ax.plot(x, yH04, '-.r', linewidth=2,
        label=r'Heyer \& Brunt 2004 $\sigma \propto R^{0.56}$')

# Bolatto+08: sigma = 0.44 * R^0.6 km/s
yB08 = 0.44 * x**0.60
ax.plot(x, yB08, ':b', linewidth=1.5,
        label=r'Bolatto+08 $\sigma \propto R^{0.60}$')

# Larson81
ax.plot(x, y, '--k', linewidth=1.5,
        label=r'Larson 1981 $\sigma \propto R^{0.38}$')

# Add literature
litpath = 'literature/'
xegc, yegc = np.loadtxt(litpath + 'ExtraGalacticGMCs.csv',  # Bolatto08
                        delimiter=',', unpack=True)
xgc, ygc = np.loadtxt(litpath + 'GalacticCenter.csv',
                      delimiter=',', unpack=True)
x64, y64 = np.loadtxt(litpath + 'M64.csv', delimiter=',', unpack=True)
# normalization ~5x higher than found in MW
xmark, ymark, ymark_err = np.loadtxt(litpath + 'SMMJ2135.txt', unpack=True)


# read in
ax.scatter(xmark, ymark, label="SMM J2135-0102", color='magenta', marker='D', s=10)
ax.scatter(x64, y64, label="M64", color='orange', marker='^', s=10)
ax.scatter(xgc, ygc, label="Heyer Galactic Center", color='b', marker='o', s=7)
ax.scatter(xegc, yegc, label="Bolatto+08: Extra-Galactic GMCs",
           color='k', marker='.', s=10)

ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].R_pc * 2.0)
        _y.append(np.sqrt(ss[ks][kkk].sigmaSq) * cm2km)

    ax.plot(_x, _y, ls='',
        markersize=10, marker='*', label='snapshot: ' + str(ks))

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Cloud Size [pc]")
ax.set_ylabel(r"$\sigma$ [km s$^{-1}$]")

# Shrink current axis by 20%
box = ax.get_position()
# ax.set_position([box.x0, box.y0 + 0.25, box.width, box.height])
ax.legend(loc="upper center", ncol=4, fontsize=9, bbox_to_anchor=(0.5, -0.18))

ax.set_xlim(1.0, 1.0e3)
ax.set_ylim(0.1, 7.e2)

plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'LarsonsLike_plot.png', bbox_inches="tight")


# --- sigma_v - gas mass SD plot ---
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].R_pc * 2.0)
        _y.append(ss[ks][kkk].massSD)

    ax.plot(_x, _y, ls='', markersize=10,
            marker='*', label='snapshot: ' + str(ks))

# # overplot Sun+18 relation from PHANGS+M51
# x = np.linspace(18.0, 2000.0, 1000.0)
# y = (x/1.e2)**0.47 + 10**0.85
# ax.plot(x, y, '--k', label="Sun+18: PHANGS+M51")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlim(1, 1.e5)
ax.set_ylim(1.0, 1.e4)

ax.set_xlabel(r"$\Sigma_{\rm gas}$ [M$_{\odot}$ pc$^{-2}$]")
ax.set_ylabel(r"$\sigma$ [km s$^{-1}$]")

ax.legend(loc="best")
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'veldisp_gasSD.png', bbox_inches="tight")


# --- Mgas - R^2 relation --
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].R_pc * 2.0)
        _y.append(ss[ks][kkk].mass_Msun)
    ax.plot(_x, _y, ls='', markersize=10,
            marker='*', label='snapshot: ' + str(ks))

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlim(1, 1.e3)
ax.set_ylim(1.0e3, 1.0e8)

ax.set_xlabel(r"$R^2$ [pc$^{2}$]")
ax.set_ylabel(r"$M_{\rm cl}$ [M$_\odot$]")

ax.legend(loc="best")
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'Mass_R2.png', bbox_inches="tight")


# --- R - Tff, color by Mcl --
from matplotlib.ticker import NullFormatter
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    _z = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].tff_Myr)
        _y.append(ss[ks][kkk].R_pc)
        _z.append(ss[ks][kkk].mass_Msun)

    cax = ax.scatter(_x, _y, s=np.array(_z)/np.array(_z).min(),
                     marker='*',
                     label='snapshot: ' + str(ks))

# ax.set_xscale("log")
ax.set_yscale("log")

ax.set_ylim(10.0, 1.e2)
ax.set_xlim(1.0, 5.0)

# ax.yaxis.set_major_locator(plt.MaxNLocator(1))
# ax.yaxis.set_minor_locator(plt.MaxNLocator(1))
ax.yaxis.set_minor_formatter(NullFormatter())

ax.set_xlabel(r"$t_{\rm ff}$ [Myr]")
ax.set_ylabel(r"$R$ [pc]")

ax.legend(loc="best")
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'R-tff-Mass.png', bbox_inches="tight")


# --- PVE diagram ---
fig = plt.figure()
ax = fig.add_subplot(111)

import astropy.constants as C
# V_0^2 = sigma^2/R = 1/3 (pi * Gamma * G * Sigma + 4 * P / Sigma)
# K cm^-3, Elmegreen+89: typical Pressure for neutral ISM ~9 x 1E3 K
# cm^-3, whereas
POverKb = [1.e4, 1.e5, 1.e6, 1.e7]
sd_cgs = np.logspace(-3, 1.0, 30)

# def V0_sq_func(pressure, sd_cgs, Gamma=3./5):
#     """

#     Gamma = 3./5 for uniform density sphere

#     """
#     print np.pi * Gamma * C.G.cgs.value * sd_cgs
#     print pressure/sd_cgs
#     import pdb; pdb.set_trace()

#     V0_sq = 1./3 * (np.pi * Gamma * C.G.cgs.value * sd_cgs + 4. * pressure/sd_cgs)
#     return V0_sq

# Plabel = ['4', '5', '6', '7']
# cm2pc = 1./3.086E+18
# for ip, Pext in enumerate(POverKb):
#     V0_sq = V0_sq_func(Pext, sd_cgs)
#     ax.plot(sd_cgs, V0_sq/(1.e5)**2 / cm2pc,
#             label=r'Log P = ' + Plabel[ip] + r"K cm$^{-2}$")

# Heyer+09 GRS data
xxx, yyy = np.loadtxt(litpath + "GRS.txt", unpack=True)
ax.scatter(xxx, yyy, marker='.', color='k', s=7, label='Heyer+09 GRS')

ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ks in iter(sorted(ss.iterkeys())):
    _x = []
    _y = []
    for kkk in ss[ks].iterkeys():
        _x.append(ss[ks][kkk].massSD *
              ss[ks][kkk].Msun2g / (ss[ks][kkk].pc2cm)**2)
        _y.append(ss[ks][kkk].sigmaSq / (1.e5)**2 / ss[ks][kkk].R_pc)

    ax.plot(_x, _y, ls='', markersize=10,
            marker='*', label='snapshot: ' + str(ks))

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlim(10**-3, 10**1.5)
ax.set_ylim(10**-1.5, 10**2.5)
# ax.yaxis.set_minor_formatter(NullFormatter())

ax.set_xlabel(r"$\Sigma_{\rm gas}$ [g cm$^{-2}$]")
ax.set_ylabel(r"$\sigma^2/R$ [km$^2$ s$^{-2}$ pc$^{-1}$]")

ax.legend(loc="best")
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'PVE-like.png', bbox_inches="tight")


# --- Cloud mass distribution, of clouds ID across all snapshots -----
fig = plt.figure()
ax = fig.add_subplot(111)

allmasses = []
for snap in ss.iterkeys():
    for snapleafs in ss[snap].iterkeys():
        allmasses.append(ss[snap][snapleafs].mass_Msun)
allmasses = np.array(allmasses)

# unbinned CDF
X2 = np.sort(allmasses)
F2 = np.array(range(len(allmasses)))/float(len(allmasses))
ax.plot(X2, F2, label='unbinned CDF', lw=2, alpha=1, zorder=2)

# binned CDF
H, X1 = np.histogram(allmasses, bins=10, normed=True)
dx = X1[1] - X1[0]
F1 = np.cumsum(H) * dx
ax.plot(X1[1:], F1, label='binned CDF, bins=10', lw=2.0, alpha=0.9, zorder=3)

H, X1 = np.histogram(allmasses, bins=50, normed=True)
dx = X1[1] - X1[0]
F1 = np.cumsum(H) * dx
ax.plot(X1[1:], F1, label='binned CDF, bins=50', lw=2.0, alpha=0.9, zorder=3)

H, X1 = np.histogram(allmasses, bins=100, normed=True)
dx = X1[1] - X1[0]
F1 = np.cumsum(H) * dx
ax.plot(X1[1:], F1, label='binned CDF, bins=100', lw=2.0, alpha=0.9, zorder=3)

H, X1 = np.histogram(allmasses, bins=200, normed=True)
dx = X1[1] - X1[0]
F1 = np.cumsum(H) * dx
ax.plot(X1[1:], F1, label='binned CDF, bins=200', lw=2.0, alpha=0.9, zorder=3)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$M_{\rm cl}$ [M$_\odot$]")
ax.set_ylabel("CDF")

ax.set_xlim(allmasses.min(), allmasses.max())

ax.legend(loc="best")
plt.tight_layout()
plt.show()
fig.savefig(leafdir_out + 'MassDistribution.png', bbox_inches="tight")


# -----------
