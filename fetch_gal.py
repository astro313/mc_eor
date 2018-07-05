"""

Get density, delta x, x_vector in subregion defined by Andrea's .csv file

last mod: 5 July 2018

"""

import pymses
from pymses.filters import PointFunctionFilter
from pymses.utils import constants as C
import os
from pymses.analysis.visualization import Camera
from numpy import array
import numpy as np
import pylab
from pymses.analysis.visualization import *
import matplotlib.pyplot as plt

ro = pymses.RamsesOutput("output", 28)

boxlen_pc = ro.info['unit_length'].express(C.pc)
finest_res = boxlen_pc / 2**ro.info['levelmax']
# 32.09690179793066 pc

amr = ro.amr_source(['rho'])


def amr2cell(ro=None, list_var=None, log_sfera=False, camera_in={}, verbose=False):
    """
    log_sfera: Boolean
        True for sphere
    """
    assert ro != None
    assert list_var != None

    from pymses.utils import regions
    from pymses.filters import RegionFilter, CellsToPoints

    amr = ro.amr_source(list_var)

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
    amr = RegionFilter(regione_sp, amr)
    # get everithing in the region
    amr = CellsToPoints(amr)

    celle = amr.flatten()
    amr = None

    return celle


center = [0.53103, 0.51031000000000004, 0.50402000000000002]
region_size = [0.0015, 0.0015]
cells_inside_camera = amr2cell(ro, list_var=['rho'], log_sfera=False,
                               camera_in={'center': center, 'region_size': region_size}, verbose=True)

dx_vector = cells_inside_camera.get_sizes()
print len(dx_vector)

loc_vector = cells_inside_camera.points
print loc_vector.shape

# check location
print 'xmin, xmax', loc_vector[:, 0].min(), loc_vector[:, 0].max()

plt.figure('1')
density_vector = cells_inside_camera['rho']
plt.hist(np.log10(density_vector))
plt.show()


plt.hist(loc_vector[:, 0])            # denser towards center of galaxy.
plt.show()

plt.hist(np.log2(1./dx_vector))       # shows the level of refinement
plt.show()


np.savez_compressed('snapshot28_center', density_vector=density_vector, loc_vector=loc_vector, dx_vector=dx_vector)