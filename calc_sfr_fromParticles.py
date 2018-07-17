'''

Calculate the globally integrated SFR for each snapshot. (using pymses, not fully functional). Always load in "level"!!

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


def calculate_age_stars(ro_in=None, dset_in=None, time_proper=True):

    if(time_proper):
        # switch depends on ramses run setup
        import cosmolopy.distance as cd
        import cosmolopy.constants as cc

        cosmo = {'omega_M_0': ro_in.info["omega_m"],
                 'omega_lambda_0': ro_in.info["omega_l"],
                 'h': ro_in.info["H0"] / 100.
                 }
        cosmo = cd.set_omega_k_0(cosmo)

        t_z0 = cd.age(0., **cosmo) / (cc.Gyr_s /
                                      1.e+3)                         # Myr
        ram2myr = ro_in.info["unit_time"].express(
            C.Myr) / ro_in.info["aexp"]**2  # Myr
        # age of the universe when the star particle was created
        star_age = t_z0 + dset_in["epoch"][:] * ram2myr
    else:
        Myr_unit_time = ro_in.info["unit_time"].express(C.Myr)
        stars_age = (ro_in.info["time"] - dset_in["epoch"][:]) * Myr_unit_time

    return stars_age


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


def particles2cell(ro=None, star=True, list_var=None, log_sfera=False, camera_in={}, verbose=False):
    """
    log_sfera: Boolean
        True for sphere
    """

    assert ro != None
    assert list_var != None

    from pymses.utils import regions
    from pymses.filters import RegionFilter, CellsToPoints

    part = ro.particle_source(list_var)

    # Filter all the particles which are initially present in the simulation
    from pymses.filters import PointFunctionFilter
    if star:
        star_filter = lambda dset: dset["epoch"] != 0.0
        part = PointFunctionFilter(star_filter, part)
    else:
        dm_filter = lambda dset: dset["epoch"] == 0.0
        part = PointFunctionFilter(dm_filter, part)

    center = camera_in['center']
    radius = camera_in['region_size'][0]

    if(log_sfera):
        regione_sp = regions.Sphere(center, radius)
    else:
        sinistra = np.copy(center) - radius
        destra = np.copy(center) + radius
        regione_sp = regions.Box((sinistra, destra))

    if(verbose):
        print 'Extracting cells'
        if(log_sfera):
            print '  getting a sphere'
            print '  center:', center
            print '  radius:', radius
        else:
            print '  getting a box'
            print '  center:', center
            print '  size  :', radius
            print '  left  :', sinistra
            print '  right :', destra

    # cut the region
    part = RegionFilter(regione_sp, part)
    celle = part.flatten()

    return celle, part


if __name__ == '__main__':

    debug = False
    import cPickle as pickle
    import matplotlib.pyplot as plt

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

        parts_inside_camera_vec, parts_inside_camera = particles2cell(ro, star=True,
                                                                      list_var=part_fields,
                                                                      camera_in={'center': center,
                                                                                 'region_size': region_size},
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

        # position automatically loadded in.
        print parts_inside_camera_vec.points.shape

        # ----------------------------------------------------------------
        # each star particle has a mass and age. Select those w/in delta_t
        lookbackTimeGyr = parts_inside_camera_vec[
            "epoch"] * ro.info['unit_time'].express(C.Gyr)

        print (13.7 - abs(lookbackTimeGyr)
               ).max(), (13.7 - abs(lookbackTimeGyr)).min()
        # # exclude IC particles
        # idx_within10Myr = parts_inside_camera["epoch"] <= 10.0

        # convert "epoch" to Myr?
        age_star_myr = calculate_age_stars(
            ro_in=ro, dset_in=parts_inside_camera["epoch"])

        # exclude IC particles
        idx_within10Myr = parts_inside_camera["epoch"] <= 10.0

        # # convert stellar mass in code unit to Msun
        # Mstar_Msun = something here

        # SFR = Mstar_Msun[idx]
