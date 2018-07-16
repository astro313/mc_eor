"""

Get density, delta x, x_vector in subregion defined by Andrea's .csv file

+ other useful fields such as:


'x-velocity'
'y-velocity'
'z-velocity'
'Density'
'H2'
'Pressure'
'Pressure_nt'
'Z'

last mod: 16 July 2018


"""

from io_module.manipulate_fetch_gal_fields import get_units

# print (__doc__)

import pymses
from pymses.sources.ramses import output
from pymses.utils import constants as C
import numpy as np
from pymses.analysis.visualization import *
import matplotlib.pyplot as plt


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
    amr = CellsToPoints(amr)

    celle = amr.flatten()
    amr = None

    return celle


def getpoints4fields(ro, outname, fields, center, region_size, log_sfera=False, debug=False):
    """

    Parameters
    ----------
    outname: str
        output file name in .npz file format (extension will be automatically added)

    fields: list of strings
        list of strings to fetch from sim. output to save (and care to resample)

    center: list of 3 floats
        center on which to extract points

    region_size: list of floats
        radius if sphere

    log_sfera: bool
        if True: reigon is a sphere, else it's a box

    debug: bool

    Returns
    -------
    None

    """

    cells_inside_camera = amr2cell(ro,
                                   list_var=fields,
                                   log_sfera=log_sfera,
                                   camera_in={'center': center,
                                              'region_size': region_size},
                                   verbose=debug)

    dx_vector  = cells_inside_camera.get_sizes() # width of the cells
    loc_vector = cells_inside_camera.points      # position of the center of the cells

    dict_unit  = get_units(ro=ro)                # unit dictionary to convert variables

    if debug:
        print "length of dx_vector: ", len(dx_vector)
        print "shape of x_vector: ", loc_vector.shape
        # check location
        print 'xmin, xmax', loc_vector[:, 0].min(), loc_vector[:, 0].max()

        plt.figure()

        for ii in fields:

            norm, unit_lab = 1., ''
            try:
                norm, unit_lab = dict_unit[ii][0], dict_unit[ii][1]
            except KeyError:
                pass

            _vector = cells_inside_camera[ii]

            to_plot = np.copy(_vector)
            to_plot = to_plot * norm

            plot_title = ii
            if(unit_lab != ''):
                plot_title = plot_title + ' / ' + unit_lab

            if ii in ['rho', 'P', 'P_nt']:
                to_plot = np.log10(to_plot)
                plot_title = 'log(' + plot_title + ')'

            plt.hist(to_plot)
            plt.title(plot_title)
            plt.show()

        # expect to be denser towards center of galaxy.
        plt.hist(loc_vector[:, 0])
        plt.title('x-vector (code unit)/postion of reigon extracted from big box')
        plt.show()

        plt.hist(np.log2(1. / dx_vector))
        plt.title('level of refinement')
        plt.show()

    param_dict = {'dx_vector': dx_vector, 'loc_vector': loc_vector}
    for i in fields:
        param_dict[i] = cells_inside_camera[i]

    if debug:
        print param_dict

    import os
    if os.path.isfile(outname):
        os.system('rm ' + outname)
    if debug:
        print 'Saving data to'
        print '  ', outname
    np.savez_compressed(outname, **param_dict)

    return None



if __name__ == "__main__":
    ro = pymses.RamsesOutput("output", 28)

    boxlen_pc = ro.info['unit_length'].express(C.pc)  # 32.09690179793066 pc
    finest_res = boxlen_pc / 2**ro.info['levelmax']

    center = [0.53103, 0.51031000000000004, 0.50402000000000002]
    region_size = [0.0015, 0.0015]

    fields = ['rho', 'vel', 'P_nt', 'P', 'H2', 'Z']

    camera_in = {'center': center,
                 'region_size': region_size}

    getpoints4fields(ro, 'snapshot28_center_fields0123456-15',
                     fields, center, region_size, log_sfera=False, debug=True)
