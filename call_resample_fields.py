'''

Last mod: 13 July 2018


'''

import os
import numpy as np
from resample_fields import resam_each_field, resam_vel_field
import pickle

folder = 'precomputed_data/'
f_camera = folder + 'camera_settings.log'
with open(f_camera,'rb') as f:
    data = pickle.load(f)

snapshotToLoad = range(16, 28)
fieldsToExtract = ['rho', 'P_nt', 'P', 'H2', 'Z']

for ssnum in snapshotToLoad:

   # saved in call_fetch_gal_fields.py
    ds = np.load('snapshot' + str(ssnum) + '_center_fields0123456-15.npz')
    outname = "snapshot" + str(ssnum) + "_center_fields0123456-15_resampled.h5"

    if os.path.exists(outname):
        os.system('rm ' + outname)

    originalSize = [data[str(ssnum)]['size']]

    #
    dx_vector = ds['dx_vector']
    loc_vector = ds['loc_vector']

    fieldsDict = {}
    for i in fieldsToExtract:
        fieldsDict[i] = ds[i]

    for kkk, vvv in fieldsDict.iteritems():
        resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, debug=False)

    # -- velocity --
    axes = {'vel_x': ds['vel'][:, 0], 'vel_y': ds['vel'][:, 1], 'vel_z': ds['vel'][:, 2]}

    for kkk, vvv in axes.iteritems():
      resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, originalSize, debug=True)


# ------