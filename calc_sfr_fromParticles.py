'''

Calculate the globally integrated SFR for each snapshot.

Last mod: 16 July 2018

'''


from io_modules.manipulate_fetch_gal_fields import get_units

# print (__doc__)

import pymses
from pymses.sources.ramses import output
from pymses.utils import constants as C
import numpy as np
from pymses.analysis.visualization import *


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

def particles2cell(ro=None, list_var=None, log_sfera=False, camera_in={}, verbose=False):

    """
    log_sfera: Boolean
        True for sphere
    """
    assert ro != None
    assert list_var != None

    from pymses.utils import regions
    from pymses.filters import RegionFilter, CellsToPoints

    part = ro.particle_source(list_var)

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
    part = None

    return celle


if __name__ == '__main__':

    debug = False
    import cPickle as pickle

    part_fields = ['vel', "id", "epoch", "mass"]

    delta_t = 10.0    # Myr

    folder = 'precomputed_data/'
    f_camera = folder + 'camera_settings.log'
    dist = 0.00075     # from output/output_00028/camera_28_[cut_13lev].csv
    far_cut_depth = 0.00075

    with open(f_camera, 'rb') as f:
        data = pickle.load(f)

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

        # Filter all the particles which are initially present in the simulation
        from pymses.filters import PointFunctionFilter
        dm_filter = lambda dset: dset["epoch"] == 0.0
        dm_parts = PointFunctionFilter(dm_filter, point_dset)

        parts_inside_camera = particles2cell(ro,
                                   list_var=part_fields,
                                   camera_in={'center': center,
                                              'region_size': region_size},
                                   verbose=debug)

        # visualize
        # map operator: mass
        scal_func = ScalarOperator(lambda dset: dset["mass"])     # simple, plot the mass

        # map processing
        mp = fft_projection.MapFFTProcessor(parts_inside_camera, ro.info)

        cam = Camera(center=camera_in['center'],
                     line_of_sight_axis=camera_in['los'],
                     region_size=[camera_in['region_size'][0], camera_in['region_size'][0]],
                     distance=dist,
                     far_cut_depth=far_cut_depth,
                     up_vector=camera_in['up_vec'],
                     map_max_size=camera_in['map_max_size'],
                     log_sensitive=True)
        mapp = mp.process(scal_func, cam, surf_qty=True)
        P.imshow(np.log10(mapp))
        P.show()

        # can't plot???


        # ----------------------------------------------------------------
        # each star particle has a mass and age. Select those w/in delta_t

        # convert "epoch" to Myr?
        something here...

        # exclude IC particles
        idx_within10Myr = parts_inside_camera["epoch"] <= 10.0


        # convert stellar mass in code unit to Msun
        Mstar_Msun = something here

        SFR = Mstar_Msun[idx]


