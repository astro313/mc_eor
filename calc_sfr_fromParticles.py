'''

Calculate the globally integrated SFR for each snapshot. Using pymses.

NOTE
----
Always load in "level" for AMR particles to work!


Last mod: 18 July 2018

'''


from io_modules.manipulate_fetch_gal_fields import get_units

# print (__doc__)

import pymses
from pymses.utils import constants as C
import numpy as np
import os
from pymses.analysis.visualization import *


if __name__ == '__main__':

    stardir = 'star_particles/'

    if not os.path.isdir(stardir):
        os.mkdir(stardir)

    from io_modules.pymses_helper import load_in_amr
    load_in_amr()

    debug = False
    import cPickle as pickle
    import matplotlib.pyplot as plt
    from io_modules.pymses_helper import particles2cell, calculate_age_stars

    part_fields = ['vel', "id", "epoch", "mass", 'level']

    delta_t = [4.0, 10.0, 20.0, 100.0]    # Myr
    assert len(delta_t) == 4      # because of the way we wrote I/O

    folder = 'precomputed_data/'
    f_camera = folder + 'camera_settings.log'
    dist = 0.00075     # from output/output_00028/camera_28_[cut_13lev].csv
    far_cut_depth = 0.00075

    with open(f_camera, 'rb') as f:
        data = pickle.load(f)

    starMs = []
    sfrs_MsunPerYr = {}

    snapshotsToLoad = range(16, 29)
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

        parts_inside_camera_vec,\
            parts_inside_camera =\
            particles2cell(ro, star=True, list_var=part_fields,
                                camera_in=
                                {
                               'center': center,
                               'region_size': region_size
                               },
                               verbose=debug)

        mass = parts_inside_camera_vec['mass'] * \
            ro.info['unit_mass'].express(C.Msun)

        plt.figure()
        plt.hist(mass, bins=200)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Stellar Mass [Msun]')
        plt.title("Mass of Star Particles in snapshot " + str(ssnum))
        plt.savefig(stardir + "starMassHist_" + str(ssnum) + ".png")
        plt.close()

        # visualize
        # map operator: mass
        scal_func = ScalarOperator(lambda dset: dset["mass"] * ro.info['unit_mass'].express(C.Msun) / (ro.info['unit_length'].express(C.pc))**2)     # simple, plot the mass

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
        plt.close()
        plt.figure()
        plt.imshow(np.log10(mapp))
        cbar = plt.colorbar()
        cbar.set_label(r"log10(Mass surface density [Msun/pc$^2$])")
        plt.savefig(stardir + "star_" + str(ssnum) + '.png')
        plt.close()

        # particles positions are automatically loadded in.
        print parts_inside_camera_vec.points.shape

        # ----------------------------------------------------------------
        # each star particle has a mass and age. Select those w/in delta_t
        lookbackTimeGyr = parts_inside_camera_vec[
            "epoch"] * ro.info['unit_time'].express(C.Gyr)

        UniverseAgeStarFormationMyr = calculate_age_stars(
            ro_in=ro, dset_in=parts_inside_camera_vec)

        # i.e,. GalaxyAgeFirstStarFormationMyr
        GalaxyAgeFirstStarFormationMyr = UniverseAgeStarFormationMyr - \
            np.min(UniverseAgeStarFormationMyr)

        # stellar mass buildup history
        ii = np.argsort(GalaxyAgeFirstStarFormationMyr)
        x = GalaxyAgeFirstStarFormationMyr[ii]
        y = mass[ii]
        y = np.cumsum(y)
        plt.figure()
        plt.plot(x, y)
        print "Total Stellar Mass buildup: {:.2f} x 10^10 [Msun]".format(y.max() / 1.e10)
        starMs.append(y.max() / 1.e10)

        plt.ylabel("Stellar Mass [Msun]")
        plt.xlabel("Age of Universe [Myr]")
        plt.savefig(stardir + "snapshot" + str(ssnum) + '_MstarBuildup.png')
        plt.close()

        _sfrs = {}
        for it in delta_t:
            idx_withinNMyr = np.abs(np.max(GalaxyAgeFirstStarFormationMyr) - GalaxyAgeFirstStarFormationMyr) <= it

            print 'SFR within', it, 'Myr'
            print '  ', np.sum(mass[idx_withinNMyr]) / 1.e+6 / it, 'Msun/yr'
            _sfrs[str(int(it))] = (np.sum(mass[idx_withinNMyr]) / 1.e+6 / it)

            # get SFR averaged over delta_t
            x = GalaxyAgeFirstStarFormationMyr[ii]
            x_binned_myr = np.arange(np.min(x), np.max(x), it)

            y, t_binned = np.histogram(
                x, weights=mass[ii], bins=int((x.max() - x.min()) / it))
            t_binned = t_binned[1:]

            y = y / (np.mean(np.diff(t_binned)) * 1.e6)

            plt.figure()
            plt.plot(t_binned, y)
            plt.ylabel("dM/dt [Msun/yr]")
            plt.xlabel("Age of Galaxy [Myr]")
            plt.title(" SFR of Main Galaxy in snapshot " + str(ssnum) + " averaged over " + str(int(it)) + " Myr")
            plt.tight_layout()
            plt.savefig(stardir + "snapshot" + str(ssnum) + '_SFR_' + str(int(it)) + 'Myr.png')
            plt.close()

        sfrs_MsunPerYr[str(ssnum)] = _sfrs

    starMS = np.round(np.array(starMs), 2)

    # dump out into a file
    outfile = stardir + 'SFRs.txt'
    File = open(outfile, 'w')
    File.write("# snapshot\t M*_built\t SFR_" + str(int(np.sort(delta_t)[0])) + '\t SFR_' + str(int(np.sort(delta_t)[1])) + '\t SFR_' + str(int(np.sort(delta_t)[2])) + '\t SFR_' + str(int(np.sort(delta_t)[3])))
    File.write("\n# M*_built in [10^10 Msun].")
    File.write("\n# SFR_n in [Msun/yr], averaged over n Myr.\n")
    for idx, ss in enumerate(iter(sorted(sfrs_MsunPerYr.iterkeys()))):
        File.write("{0:s}\t{1:.2f}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\n".format(ss, starMs[idx], sfrs_MsunPerYr[ss][str(int(np.sort(delta_t)[0]))], sfrs_MsunPerYr[ss][str(int(np.sort(delta_t)[1]))], sfrs_MsunPerYr[ss][str(int(np.sort(delta_t)[2]))], sfrs_MsunPerYr[ss][str(int(np.sort(delta_t)[3]))]))
    File.close()

print 'Stellar mass buildup for each snapshot: '
print ' ', np.round(np.array(starMs), 2), 'x 10^10 Msun/yr'

# ----------
