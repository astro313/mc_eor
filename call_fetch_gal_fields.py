'''

fetch AMR and star fields for snapshots 16-28

Last mod: 23 July 2018

'''


import os
import pickle
from io_modules.fetch_gal_fields import getpoints4fields, get_stars_for_fields
import pymses
from pymses.utils import constants as C


folder = 'precomputed_data/'
f_camera = folder + 'camera_settings.log'

istart = 28
istop  = 29

snapshotToLoad = range(istart,istop)

fieldsToLoad = ['rho', 'vel', 'P_nt', 'P', 'H2', 'Z']
fields_stars = ['id', 'vel', 'epoch', 'mass', 'level']

with open(f_camera, 'rb') as f:
    data = pickle.load(f)

for ssnum in snapshotsToLoad:
    ro = pymses.RamsesOutput("output", ssnum)

    boxlen_pc = ro.info['unit_length'].express(C.pc)  # 32.09690179793066 pc
    finest_res = boxlen_pc / 2**ro.info['levelmax']

    center = data[str(ssnum)]['center_init']
    region_size = [data[str(ssnum)]['size']]

    camera_in = {'center': center,
                 'region_size': region_size}
    # # AMR
    getpoints4fields(ro, 'snapshot' + str(ssnum) + '_center_fields0123456-15',
                     fieldsToLoad, center, region_size,
                     log_sfera=False, debug=False)

    # star particles
    get_stars_for_fields(ro, 'snapshot' + str(ssnum) + '_center_stars',
                     fields_stars, center, region_size,
                     log_sfera=False, debug=False, star=True)
