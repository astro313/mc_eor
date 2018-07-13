'''

Last mod: 12 July 2018


'''

import os
import numpy as np
from resample_fields import resam_each_field, resam_vel_field

snapshotToLoad = range(16, 28)
fieldsToExtract = ['rho', 'P_nt', 'P', 'H2', 'Z']

for ssnum in snapshotToLoad:

   # saved in fetch_gal_fields_16-27.py
    ds = np.load('snapshot' + str(ssnum) + '_center_fields0123456-15.npz')
    outname = "snapshot" + str(ssnum) + "_center_fields0123456-15_resampled.h5"

    import os
    if os.path.exists(outname):
        os.system('rm ' + outname)
    #
    dx_vector = ds['dx_vector']
    loc_vector = ds['loc_vector']

    fieldsDict = {}
    for i in fieldsToExtract:
        fieldsDict[i] = ds[i]

    for kkk, vvv in fieldsDict.iteritems():
        resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, debug=False)

    # -- velocity --
    vel = ds['vel']
    axes = {'x': vel[:, 0], 'y': vel[:, 1], 'z': vel[:, 2]}

    for kkk, vvv in axes.iteritems():
        resam_vel_field(dx_vector, loc_vector, vvv, kkk, outname, debug=False)


# ------