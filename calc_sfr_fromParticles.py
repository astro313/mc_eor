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
from plot_modules.plot_cloud_prop import setup_plot
from matplotlib import cm
cmap     = cm.get_cmap('magma')


if __name__ == '__main__':

    setup_plot()
    debug   = False
    verbose = True
    stardir = 'star_particles/'

    if not os.path.isdir(stardir):
        os.mkdir(stardir)

    from io_modules.pymses_helper import load_in_amr
    load_in_amr()

    import cPickle as pickle
    import matplotlib.pyplot as plt
    from io_modules.pymses_helper import particles2cell, calculate_age_stars

    part_fields = ["id", "epoch", "mass", 'level']

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

        center      = data[str(ssnum)]['center_init']
        region_size = [data[str(ssnum)]['size']]
        los         = data[str(ssnum)]['los_vec']
        up          = data[str(ssnum)]['up_vec']
        mms         = data[str(ssnum)]['mms']

        if debug:
          region_size = [0.1,0.1]

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

        if debug:
          xx = parts_inside_camera_vec.points[:,0]
          xx = xx*ro.info['unit_length'].express(C.kpc)
          print np.max(xx),np.min(xx)
          print np.max(xx)-np.min(xx)
          print np.shape(xx)
        mass = parts_inside_camera_vec['mass'] * \
            ro.info['unit_mass'].express(C.Msun)

        name_out = stardir + "starMassHist_" + str(ssnum) + ".png"
        if verbose:
          print 'Figure to '
          print name_out
        plt.figure()
        plt.hist(mass, bins=200)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Stellar Mass [Msun]')
        plt.title("Mass of Star Particles in snapshot " + str(ssnum))
        plt.savefig(name_out)
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
        plt.close('all')
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        #
        z1,z2 = 0.0, np.log10(np.max(mapp))

        size_region_kpc = camera_in['region_size'][0] * ro.info['unit_length'].express(C.kpc)
        bu = ax.imshow(np.log10(mapp),extent = [-size_region_kpc/2,size_region_kpc/2,-size_region_kpc/2,+size_region_kpc/2]
                       ,vmin = z1, vmax = z2,
                       cmap=cmap
                      )
        plt.ylabel("kpc", fontsize=18)
        plt.xlabel("kpc", fontsize=18)
        plt.tick_params(which='minor', length=4)
        #
        cbar = fig.colorbar(bu)
        cbar.set_label(r"$\log \Sigma_{\star}$" + r"[M$_{\odot}\,{\rm pc}^{-2}]$", fontsize=18)
        cbar.set_clim(z1,z2)
        plt.tight_layout()
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
        name_out = stardir + "snapshot" + str(ssnum) + '_MstarBuildup.png'
        if verbose:
          print 'Stellar mass to'
          print '  ',name_out
        ii = np.argsort(GalaxyAgeFirstStarFormationMyr)
        x = np.copy(GalaxyAgeFirstStarFormationMyr[ii])
        y = np.copy(mass[ii])
        y = np.cumsum(y)
        plt.figure()
        plt.plot(x, y)
        if verbose:
          print "Total Stellar Mass buildup: {:.2f} x 10^10 [Msun]".format(y.max() / 1.e10)
        starMs.append(y.max() / 1.e10)
        plt.ylabel("Stellar Mass [Msun]")
        plt.xlabel("Age of Universe [Myr]")
        plt.savefig(name_out)
        plt.close()

        _sfrs = {}
        for it in delta_t:
            idx_withinNMyr = np.abs(np.max(GalaxyAgeFirstStarFormationMyr) - GalaxyAgeFirstStarFormationMyr) <= it

            if verbose:
              print 'SFR within', it, 'Myr'
              print '  ', np.sum(mass[idx_withinNMyr]) / 1.e+6 / it, 'Msun/yr'
            _sfrs[str(int(it))] = (np.sum(mass[idx_withinNMyr]) / 1.e+6 / it)

            # get SFR averaged over delta_t
            if(1):
              x = GalaxyAgeFirstStarFormationMyr[ii]
              x_binned_myr = np.arange(np.min(x), np.max(x), it)
              y, t_binned = np.histogram(
                  x, weights=mass[ii], bins=int((x.max() - x.min()) / it))
              t_binned = t_binned[1:]
              y = y / (np.mean(np.diff(t_binned)) * 1.e6)
            else:
              # equivalent, only for testing
              t_binned = np.arange(np.min(GalaxyAgeFirstStarFormationMyr),np.max(GalaxyAgeFirstStarFormationMyr),it)
              print np.min(GalaxyAgeFirstStarFormationMyr),np.max(GalaxyAgeFirstStarFormationMyr)
              print t_binned
              bin_in_x = np.digitize(GalaxyAgeFirstStarFormationMyr, t_binned[:-1])
              y        = np.bincount(bin_in_x,weights = mass) / ( it * 1.e6)

            out_name = stardir + "snapshot" + str(ssnum) + '_SFR_' + str(int(it)) + 'Myr.png'
            if verbose:
              print '  SFR averaged over',it,'Myr to'
              print '    ',name_out
            plt.figure()
            plt.plot(t_binned, y)
            plt.ylabel("SFR " + r"[$M_{\odot}\,yr^{-1}$]", fontsize=18)
            plt.xlabel("Age of Galaxy [Myr]", fontsize=18)
            plt.title(" SFR of Main Galaxy in snapshot " + str(ssnum) + " averaged over " + str(int(it)) + " Myr")
            plt.tight_layout()
            plt.savefig(out_name)
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
