'''

Last mod: 23 July 2018


'''

import os
import numpy as np
from gridding_modules.resample_fields import resam_each_field
from gridding_modules.locate_parts import resample_star_mass_age
import pickle
import h5py
import matplotlib.pyplot as plt


plotdir = "resampled/"
if not os.path.isdir(plotdir):
    os.mkdir(plotdir)

folder = '/mnt/home/daisyleung/mc_eor/precomputed_data/'
f_camera = folder + 'camera_settings.log'
with open(f_camera, 'rb') as f:
    data = pickle.load(f)

snapshotToLoad = range(16, 29)
fieldsToExtract = ['rho', 'P_nt', 'P', 'H2', 'Z']

# for ssnum in snapshotToLoad:
# 
#     # saved in call_fetch_gal_fields.py
#     ds = np.load('snapshot' + str(ssnum) + '_center_fields0123456-15.npz')
#     outname = "snapshot" + str(ssnum) + "_center_fields0123456-15_resampled.h5"
# 
#     if os.path.exists(outname):
#         os.system('rm ' + outname)
# 
#     originalSize = data[str(ssnum)]['size']
# 
#     #
#     dx_vector = ds['dx_vector']
#     loc_vector = ds['loc_vector']
# 
#     fieldsDict = {}
#     for i in fieldsToExtract:
#         fieldsDict[i] = ds[i]
# 
#     for kkk, vvv in fieldsDict.iteritems():
#         N = resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, originalSize, debug=False)
# 
#     # -- velocity --
#     axes = {'vel_x': ds['vel'][:, 0], 'vel_y': ds['vel'][:, 1], 'vel_z': ds['vel'][:, 2]}
# 
#     for kkk, vvv in axes.iteritems():
#         resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, originalSize, debug=False)
# 
#     # ------ sanity check -- still in code units..
#     f = h5py.File(outname, "r")
#     # careful sometimes i used "density", see
#     # resample_fields.py to make sure
#     H2density = f["rho"].value * f["H2"].value
#     plt.figure()
#     plt.imshow(np.log10(H2density[:, :, 2**4:].sum(axis=0)))
#     plt.tight_layout()
#     plt.savefig(plotdir + str(ssnum) + "_h2density.png")
# 
#     plt.figure()
#     plt.hist(np.log10(H2density[H2density > 1.e-23].flatten()), bins=100)
#     plt.tight_layout()
#     plt.savefig(plotdir + str(ssnum) + "_h2densityHist.png")
# 

# star particles: to merge w/ above
for ssnum in snapshotToLoad:   # to merge w/ above

    # saved in call_fetch_gal_fields.py
    ds = np.load('snapshot' + str(ssnum) + '_center_stars.npz')
    outname = "snapshot" + str(ssnum) + "_center_stars_resampled.h5"

    if os.path.exists(outname):
        os.system('rm ' + outname)

    center = data[str(ssnum)]['center_init']   # to merge
    originalSize = data[str(ssnum)]['size']   # to merge

    camera = {'center': center,   # to merge
              'region_size': originalSize}
    loc_vec = ds['loc_vector']
    level_vec = ds['level']
    epoch_vec = ds['epoch']
    id_vec = ds['id']
    mass_vec = ds['mass']

    # hard-coded for now because don't want to re-run
    amrfile = 'snapshot' + str(ssnum) + '_center_fields0123456-15_resampled.h5'
    f = h5py.File(amrfile, "r")
    N = f['rho'].shape[0]

    # resample star particle mass and compute mass-weighted avg epoch correspondingly                        
    resample_star_mass_age(loc_vec, level_vec, epoch_vec, mass_vec, outname, camera, Nrefined=N, debug=True)

    # ------ sanity check -- still in code units..
    f = h5py.File(outname, "r")
    mass = f['mass'].value
    plt.figure()
    plt.imshow(np.log10(mass[:, :, :].sum(axis=0)))
    plt.tight_layout()
    plt.savefig(plotdir + str(ssnum) + "_massCollapsed.png")

    plt.figure()
    plt.hist(np.log10(mass[mass > 1.e-23].flatten()), bins=100)
    plt.tight_layout()
    plt.savefig(plotdir + str(ssnum) + "_massCollapsedHist.png")
