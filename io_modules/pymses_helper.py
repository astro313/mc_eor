import numpy as np
import pymses
from pymses.utils import constants as C


def particles2cell(ro=None, star=True, list_var=None, log_sfera=False, camera_in={}, verbose=False):
    """
    log_sfera: Boolean
        True for sphere
    """

    assert ro is not None
    assert list_var is not None

    from pymses.utils import regions
    from pymses.filters import RegionFilter

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


def calculate_age_stars(ro_in=None, dset_in=None, time_proper=True):

    """

    Return
    ------
    starsFormedatUniverseAge:
        age of the universe when the star particle was created in Myr

    """

    if(time_proper):
        # switch depends on ramses run setup
        import cosmolopy.distance as cd
        import cosmolopy.constants as cc

        cosmo = {'omega_M_0': ro_in.info["omega_m"],
                 'omega_lambda_0': ro_in.info["omega_l"],
                 'h': ro_in.info["H0"] / 100.
                 }
        cosmo = cd.set_omega_k_0(cosmo)

        t_z0 = cd.age(0., **cosmo) / (cc.Gyr_s / 1.e+3)    # Myr
        ram2myr = ro_in.info["unit_time"].express(
            C.Myr) / ro_in.info["aexp"]**2

        starsFormedatUniverseAge = t_z0 + dset_in["epoch"][:] * ram2myr
    else:
        Myr_unit_time = ro_in.info["unit_time"].express(C.Myr)
        starsFormedatUniverseAge = (ro_in.info["time"] - dset_in["epoch"][:]) * Myr_unit_time

    return starsFormedatUniverseAge
