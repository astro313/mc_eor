'''

Calculate the globally integrated SFR for each snapshot. (using yt)

Last mod: 16 July 2018

'''

import yt
from io_modules.manipulate_fetch_gal_fields import get_units
import numpy as np

from yt.analysis_modules.star_analysis.api import StarFormationRate
import matplotlib.pyplot as plt


def young_stars(pfilter, data, delta_t_Myr):
    age = data.ds.current_time - data[pfilter.filtered_type, "creation_time"]
    filter = np.logical_and(age.in_units('Myr') <= delta_t_Myr, age >= 0)
    return filter


if __name__ == '__main__':

    debug = False
    import cPickle as pickle

    part_fields = ['vel', "id", "epoch", "mass"]

    delta_t_Myr = 10.0    # Myr

    folder = 'precomputed_data/'
    f_camera = folder + 'camera_settings.log'
    dist = 0.00075     # from output/output_00028/camera_28_[cut_13lev].csv
    far_cut_depth = 0.00075

    with open(f_camera, 'rb') as f:
        data = pickle.load(f)

    snapshotsToLoad = range(16, 29)
    snapshotsToLoad = range(28, 29)
    for ssnum in snapshotsToLoad:
        ds = yt.load("output/output_000" + str(ssnum) + "/info_000" + \
                    str(ssnum) + ".txt")
        print ds.derived_field_list
        print ds.field_list

        # [('all', 'particle_age'), ('all', 'particle_identifier'), ('all', 'particle_mass'), ('all', 'particle_metallicity'), ('all', 'particle_position_x'), ('all', 'particle_position_y'), ('all', 'particle_position_z'), ('all', 'particle_refinement_level'), ('all', 'particle_velocity_x'), ('all', 'particle_velocity_y'), ('all', 'particle_velocity_z'), ('io', 'particle_age'), ('io', 'particle_identifier'), ('io', 'particle_mass'), ('io', 'particle_metallicity'), ('io', 'particle_position_x'), ('io', 'particle_position_y'), ('io', 'particle_position_z'), ('io', 'particle_refinement_level'), ('io', 'particle_velocity_x'), ('io', 'particle_velocity_y'), ('io', 'particle_velocity_z'), ('ramses', 'Density'), ('ramses', 'Metallicity'), ('ramses', 'Pressure'), ('ramses', 'var12'), ('ramses', 'var13'), ('ramses', 'var14'), ('ramses', 'var15'), ('ramses', 'var16'), ('ramses', 'x-Bfield-left'), ('ramses', 'x-Bfield-right'), ('ramses', 'x-velocity'), ('ramses', 'y-Bfield-left'), ('ramses', 'y-Bfield-right'), ('ramses', 'y-velocity'), ('ramses', 'z-Bfield-left'), ('ramses', 'z-Bfield-right'), ('ramses', 'z-velocity')]

        # is 'io' refering to "stars"??


        # for field in ds.derived_field_list:
        #     if field[0] == 'Stars':  # 'io'
        #         print (field)

        for field in ds.field_list:
            if field[0] == 'io':
                print (field)

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

        sp = ds.sphere(camera_in['center'], camera_in['region_size'][0])
# r        sfr = StarFormationRate(ds, data_source=sp)

        mass = sp[("io", "particle_mass")].in_units('Msun')
        lookBack_Gyr = sp[("io", "age")].in_units('Gyr')   # is this age of stellar particle in lookback time?

        from yt.utilities.cosmology import Cosmology
        co = Cosmology()
        Gyr2s = 3.15576E+16

        ageUniverse_z0_Gyr = co.hubble_time(0.).in_units("Gyr")
        ageUniverseAtz_Gyr = ageUniverse_z0_Gyr - abs(lookBack_Gyr)
        print ageUniverseAtz_Gyr.max()
        print(co.z_from_t(ageUniverseAtz_Gyr*Gyr2s)).max()
        print(co.z_from_t(1.0*Gyr2s))   # z=6 ~1 Gyr

        # assert age

        # yt.add_particle_filter("young_stars", function=young_stars, filtered_type='Stars', requires=["creation_time"])

        # ds.add_particle_filter('young_stars')

        # # mass = sp[("stars", "particle_mass")].in_units('Msun')
        # # age = sp[("stars", "age")].in_units('Myr')
        # # ct = sp[("stars", "creation_time")].in_units('Myr')

        # # Pick out only stars created w/in past 10 Myr
        # threshold = ds.quan(delta_t_Myr, "Myr")
        # mass_young = mass[age <= threshold]
        # ct_young = ct[age <= threshold]


        # sfr = StarFormationRate(ds, star_mass=mass_young,
        #                         star_creation_time=ct_young,
        #                         volume=sp.volume())

        # sfr.write_out(name="SFR_snapshot" + str(ssnum) + ".out")

        # plt.plot(sfr.lookback_time, sfr.Msol_yr)
        # plt.xlabel('Lookback Time [yr]')
        # plt.ylabel(r'SFR [M_\odot yr$^{-1}$]')
        # plt.title("SFR of main galaxy in snapshot #" + str(ssnum))
        # plt.show()

        # print "SFR over past {:.2f} Myr= {:.2f}".format(delta_t, mass_young.sum()/delta_t)


# =========
