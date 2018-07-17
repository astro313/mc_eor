'''

Calculate the globally integrated SFR for each snapshot. Using pymses.

NOTE
----
Always load in "level" for AMR particles to work!


Last mod: 17 July 2018

'''


from io_modules.manipulate_fetch_gal_fields import get_units

# print (__doc__)

import pymses
from pymses.sources.ramses import output
from pymses.utils import constants as C
import numpy as np
import os
from pymses.analysis.visualization import *


stardir = 'star_particles/'

if not os.path.isdir(stardir):
    os.mkdir(stardir)


pymses.RamsesOutput.amr_field_descrs_by_file = \
    {"2D": {"hydro": [output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]),
                      output.Vector("Bl", [4, 5, 6]),
                      output.Vector("Br", [7, 8, 9]),
                      output.Scalar("P", 10),
                      output.Scalar("Z", 11)],
            "grav": [output.Vector("g", [0, 1, 2])]},
     "3D": {"hydro": [output.Scalar("rho", 0), output.Vector("vel", [1, 2, 3]),
                      output.Scalar("P_nt", 4), output.Scalar("P", 5),
                      output.Scalar("Z", 6),
                      # # note field 7 is skipped here because it's just flags for structure of the AMR, and pymses is not picking about skipping fields
                      # output.Scalar("H", 8),
                      # output.Scalar("E", 9),
                      # output.Scalar("H+", 10),
                      # output.Scalar("HE", 11),
                      # output.Scalar("HE+", 12),
                      # output.Scalar("HE++", 13),
                      # output.Scalar("H-", 14),
                      output.Scalar("H2", 15)
                      #                      ,output.Scalar("H2+", 16)
                      ],
            "grav": [output.Vector("g", [0, 1, 2])]}}


if __name__ == '__main__':

    debug = False
    import cPickle as pickle
    import matplotlib.pyplot as plt
    from io_modules.pymses_helper import particles2cell, calculate_age_stars

    part_fields = ['vel', "id", "epoch", "mass", 'level']

    delta_t = 10.0    # Myr

    folder = 'precomputed_data/'
    f_camera = folder + 'camera_settings.log'
    dist = 0.00075     # from output/output_00028/camera_28_[cut_13lev].csv
    far_cut_depth = 0.00075

    with open(f_camera, 'rb') as f:
        data = pickle.load(f)

    snapshotsToLoad = range(28, 29)
    for ssnum in snapshotsToLoad:
        ro = pymses.RamsesOutput("output", ssnum)

        center = data[str(ssnum)]['center_init']
        region_size = [data[str(ssnum)]['size']]
        los = data[str(ssnum)]['los_vec']
        up = data[str(ssnum)]['up_vec']
        mms = data[str(ssnum)]['mms']

        camera_in = {'center': center,
                     'region_size': region_size,
                     'los': los,
                     'up_vec': up,
                     'map_max_size': mms}

        parts_inside_camera_vec, parts_inside_camera = particles2cell(ro,
                                        star=True, list_var=part_fields, camera_in={'center': center,
                                                       'region_size':region_size},
                                            verbose=debug)

        mass = parts_inside_camera_vec['mass'] * \
            ro.info['unit_mass'].express(C.Msun)
        plt.hist(mass)
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(stardir + "starMassHist_" + str(ssnum) + ".png")

        # visualize
        # map operator: mass
        scal_func = ScalarOperator(lambda dset: dset["mass"] * ro.info['unit_mass'].express(
            C.Msun) / (ro.info['unit_length'].express(C.pc))**2)     # simple, plot the mass

        # map processing
        mp = fft_projection.MapFFTProcessor(parts_inside_camera, ro.info)

        cam = Camera(center=camera_in['center'],
                     line_of_sight_axis=camera_in['los'],
                     region_size=[camera_in['region_size']
                                  [0], camera_in['region_size'][0]],
                     distance=dist,
                     far_cut_depth=far_cut_depth,
                     up_vector=camera_in['up_vec'],
                     map_max_size=camera_in['map_max_size'],
                     log_sensitive=True)
        mapp = mp.process(scal_func, cam, surf_qty=True)
        plt.imshow(np.log10(mapp))
        plt.savefig(stardir + "star_" + str(ssnum) + '.png')

        # particles positions are automatically loadded in.
        print parts_inside_camera_vec.points.shape

        # ----------------------------------------------------------------
        # each star particle has a mass and age. Select those w/in delta_t
        lookbackTimeGyr = parts_inside_camera_vec[
            "epoch"] * ro.info['unit_time'].express(C.Gyr)

        age_star_myr = calculate_age_stars(
            ro_in=ro, dset_in=parts_inside_camera_vec)

        idx_within10Myr = np.abs(np.max(age_star_myr)-age_star_myr) <= delta_t
        print 'SFR within', delta_t, 'Myr'
        print '  ', np.sum(mass[idx_within10Myr])/1.e+6/delta_t,'Msun/yr'

        # "SFH"
        ii = np.argsort(age_star_myr)
        x = age_star_myr[ii]
        y = mass[ii]
        y = np.cumsum(y)
        plt.figure()
        plt.plot(x, y)
        plt.savefig('12345.png')

        # # convert stellar mass in code unit to Msun
        # Mstar_Msun = something here

        # SFR = Mstar_Msun[idx]
