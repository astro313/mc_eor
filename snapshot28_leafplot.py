'''

plot physical properties of clouds.

last mod: 10 July 2018

NOTE
----


'''

import matplotlib
matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                            , 'font.size': 22      # points
                            , 'legend.fontsize': 16      # points
                            , 'lines.linewidth': 2       # points
                            , 'axes.linewidth': 2       # points
                            , 'text.usetex': True    # Use LaTeX to layout text
                            , 'font.family': "serif"  # Use serifed fonts
                            , 'xtick.major.size': 6     # length, points
                            , 'xtick.major.width': 2     # points
                            , 'xtick.minor.size': 3     # length, points
                            , 'xtick.minor.width': 1     # points
                            , 'ytick.major.size': 6     # length, points
                            , 'ytick.major.width': 2     # points
                            , 'ytick.minor.size': 3     # length, points

                            , 'ytick.minor.width': 1     # points
                            , 'font.serif': ("Times", "Palatino", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                            })
import matplotlib.pyplot as plt
import cPickle as pickle
from snapshot28_leafprop import Cloud
import numpy as np
import os

# calculate properties for all clouds in this snapshot, then plot.

# load back in fields of each leaf
snapshot_num = 28
leafdir = 'leaf_fields_' + str(snapshot_num) + '/'

# defined as ... '{0:.2f}'.format(args.ncut) + '_' + str(args.step) + '_'
# + str(args.Nmin) from snapshot28_leafprop.py
fname = "0.03_5_3.p"
leaf_fields = pickle.load(open(leafdir + fname, "rb"))

# hard-code to get dx for now...
# saved in fetch_gal_fields.py
from fetch_gal_fields import get_units
import pymses
ro = pymses.RamsesOutput("output", 28)
ds = np.load('snapshot28_center_fields0123456-15.npz')
dx_vector = ds['dx_vector']

originalLevels = np.log2(1. / np.unique(dx_vector))
highestRes = 2.**(originalLevels.max() * -1)
# after finely sample unigrid
dx_pc = highestRes * get_units(ro=ro)['dx'][0]


Cloud_dict = {}
for kkk in leaf_fields.iterkeys():
    Cloud_dict[kkk] = Cloud(dx_pc, leaf_fields[kkk], int(kkk))
    print Cloud_dict[kkk]

# plotting
fig = plt.figure()
ax = fig.add_subplot(111)
for kkk in leaf_fields.iterkeys():
    MMj = Cloud_dict[kkk].mass_Msun / Cloud_dict[kkk].M_jeans
    ax.scatter(Cloud_dict[kkk].mass_Msun, MMj, s=70, alpha=0.7, c='r')

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Cloud Mass [Msun]", fontsize=20)
ax.set_ylabel(r"$M_{\rm cl} / $M$_J$", fontsize=20)
plt.tight_layout()
plt.show()
fig.savefig(leafdir + 'my_cloud_sample.png', bbox_inches="tight")



fig = plt.figure()

ax = fig.add_subplot(111)
for kkk in leaf_fields.iterkeys():

    ax.scatter(Cloud_dict[kkk].mass_Msun, Cloud_dict[kkk].alpha, s=70,
               alpha=0.7, c='g')

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Cloud Mass [Msun]")
ax.set_ylabel(r"$\alpha_{\rm vir}$")

plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)

for kkk in leaf_fields.iterkeys():
    ax.scatter(Cloud_dict[kkk].mass_Msun, np.array(Cloud_dict[kkk].SFR) / 1.0e6, label="SFR KMM06", s=80, alpha=0.7, c='r', marker="*")
    ax.scatter(Cloud_dict[kkk].mass_Msun, np.array(Cloud_dict[kkk].SFR_JML) / 1.0e6 , label="SFR JML", s=80, alpha=0.7, c='r', marker="*")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Cloud Mass [Msun]")
ax.set_ylabel(r"Star Formation Rate [Msun yr$^{-1}$]")

ax.legend(loc=2)

# ax.set_xlim(1.0e3, 1.0e7)
# ax.set_ylim(1.0e-5, 1.0e1)

plt.show()
fig.savefig(leafdir + "SFR_CloudMass.png", bbox_inches="tight")


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

ax.plot(x, y, '-b', linewidth=2, label="Kennicut")
ax.plot(Heinerman_SigmaGas, Heinerman_SigmaSFR, 'bo',
        linewidth=2, alpha=0.7, label="Heinerman 2010")

for kkk in leaf_fields.iterkeys():

    ax.scatter(np.array(Cloud_dict[kkk].massSD),
               np.array(Cloud_dict[kkk].SFR) / 1.0e6 / \
               (np.pi * Cloud_dict[kkk].R_pc * pc2kpc)**2,
               label="This Work", s=80, alpha=0.7, c='r', marker='*')


ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$\Sigma_{\rm gas}$ [Msun pc$^2$]", fontsize=20)
ax.set_ylabel(r"$\Sigma_{\rm SFR}$ [Msun yr$^{-1}$ kpc$^2$]", fontsize=20)

ax.legend(loc=4)
# ax.set_xlim(5.0e0, 1.0e3)
plt.show()
fig.savefig(leafdir + 'SurfSFR_SurfGas.png', bbox_inches="tight")


# Larson's relations
fig = plt.figure()
ax = fig.add_subplot(111)

cm2km = 1.e-5

x = np.logspace(1, 3, 10)
#y = 1.3e-2*x**0.38
y = 1.01 * x**0.38

ax.plot(x, y, '--k', linewidth=2, label='Larson $\sigma \propto L^{0.38}$')

for kkk in leaf_fields.iterkeys():
    ax.scatter(Cloud_dict[kkk].R_pc * 2.0,
               np.sqrt(Cloud_dict[kkk].sigmaSq) * cm2km,
               c='k', marker='o', alpha=0.6, s=30)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Cloud Size [pc]", fontsize=20)
ax.set_ylabel(r"velocity dispersion $\sigma$ [km/s]", fontsize=20)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

ax.legend(loc=1, fontsize=20)
# ax.set_xlim(1.0e1, 1.0e3)

plt.show()
fig.savefig(leafdir + 'LarsonsLike_plot.png', bbox_inches="tight")

    # -----------
