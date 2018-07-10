'''

resample AMR points to finest resolution, do so for other useful fields in addition to "density"

which ones do we need?

'x-velocity'
'y-velocity'
'z-velocity'
'Density'
'H2'
'Pressure'
'Pressure_nt'

'cell_volume'    << from yt only?
'cell_mass'   << from yt only?
'averaged_density'  << from yt only?
'temperature'

last mod: 9 July 2018


'''

import numpy as np

# saved in fetch_gal_fields.py
ds = np.load('snapshot28_center_fields0123456-15.npz')
outname = "snapshot28_center_fields0123456-15_resampled.h5"

#
dx_vector = ds['dx_vector']
loc_vector = ds['loc_vector']

fields = ['rho', 'P_nt', 'P', 'H2', 'Z']
# extract 'vel' separately below ..

import os
if os.path.exists(outname):
    os.system('rm ' + outname)


def resam_each_field(dx_vector, loc_vector, field_vector, fieldname, outname, originalSize=0.0015, debug=True):

    """One liner description

    Parameters
    ----------
    dx_vector: array

    loc_vector: array

    field_vector: array
        subregion

    fieldname: str
        just a str (for creating dataset in h5py)

    outname: str
        output .h5 filename

    originalSize: float
        region of which we extracted this subset of fields from the original box, in code unit

    Returns
    -------

    """

    originalLevels = np.log2(1. / np.unique(dx_vector))
    numLevel = len(originalLevels)   # or numLevel = len(np.unique(dx_vector))
    levels = originalLevels               # let's just keep as
    originalLevels_vector = np.log2(1. / dx_vector)
    (originalLevels_vector == np.sort(originalLevels_vector)).all() == True

    levels_ii = originalLevels_vector
    levels_ii = np.array(levels_ii, dtype=int)
    levels = np.array(levels, dtype=int)

    d = {}
    for LL in levels:
        # find indices where they correspond to each of the levels
        mask = levels_ii == LL
        d[str(LL)] = mask
        print LL, np.max(levels_ii[mask]), np.min(levels_ii[mask])

    # check originalLevels_vector[ii] is properly saved to the right dict key
    ilowestLevel = d[str(int(levels.min()))]

    # coarset level
    Llow = str(int(levels.min()))
    pos = loc_vector[d[Llow], :]
    II = field_vector[d[Llow]]

    ilev = levels_ii[d[Llow]]
    assert ilev.min() == ilev.max() == int(Llow)

    print field_vector
    print II

    highestRes = 2.**(levels.max() * -1)
    N = originalSize / highestRes
    print N

    # make sure it's even number
    # N = 196
    N = int(N)
    if N % 2:  # odd number
        N += 1
    assert N % 2 == 0

    Ninit = N / len(levels) / 2
    field_cube = np.zeros((Ninit, Ninit, Ninit))

    imatrix = np.ones((2, 2, 2))

    for lll in np.sort(levels):
        if debug:
            print lll

        field_cube = np.kron(field_cube, imatrix)

        pos = loc_vector[d[str(lll)], :]

        xx = (pos[:, 0] - loc_vector[:, 0].min()) / \
            (loc_vector[:, 0].max() - loc_vector[:, 0].min())
        yy = (pos[:, 1] - loc_vector[:, 1].min()) / \
            (loc_vector[:, 1].max() - loc_vector[:, 1].min())
        zz = (pos[:, 2] - loc_vector[:, 2].min()) / \
            (loc_vector[:, 2].max() - loc_vector[:, 2].min())

        if debug:
            print '  x:', xx.min(), xx.max()
            print '  y:', yy.min(), yy.max()
            print '  z:', zz.min(), zz.max()
            #import pdb; pdb.set_trace()

        xx = xx * field_cube.shape[0]
        yy = yy * field_cube.shape[1]
        zz = zz * field_cube.shape[2]


        xpos = np.array(xx, dtype=int)
        ypos = np.array(yy, dtype=int)
        zpos = np.array(zz, dtype=int)

        xpos[xpos == field_cube.shape[0]] = field_cube.shape[0] - 1
        ypos[ypos == field_cube.shape[1]] = field_cube.shape[1] - 1
        zpos[zpos == field_cube.shape[2]] = field_cube.shape[2] - 1
        field_cube[xpos, ypos, zpos] = field_vector[d[str(lll)]]

        if debug:
            print '  field    :', field_cube.min(), field_cube.max()
            print '  field log:', np.log10(field_cube.min()), np.log10(field_cube.max())

    if debug:
        size = field_cube.shape[0]
        print "1d length of final resampled cube: ", size

        import matplotlib.pyplot as plt
        # 2**4: to chop off spurious edges
        if fieldname == 'rho' or fieldname == 'P_nt' or fieldname == 'P':
            plt.imshow(np.log10(field_cube[:, :, 2**4:].sum(axis=0)))
            plt.title(fieldname)
            plt.show()

        if fieldname != 'vel':
            plt.imshow(field_cube[:, :, 2**4:].sum(axis=0))
            plt.title(fieldname)
            plt.show()
        else:
            pass

    import h5py
    # if not exist, then create the .h5 file, else append mode
    if os.path.exists(outname):
        mode = "a"
    else:
        mode = "w"

    f = h5py.File(outname, mode)
    # create a dataset at the root
    f.create_dataset("/" + fieldname, data=field_cube)
    f.close()

    if debug:
        # read the h5 file back to make sure appended field
        f = h5py.File(outname, "r")
        print f.keys()
        f.close()

    return None


# --- repeat for the velocity vector ---
def resam_vel_field(dx_vector, loc_vector, field_vector, axisname, outname,
                    originalSize=0.0015, debug=True):

    """One liner description

    Parameters
    ----------
    dx_vector: array

    loc_vector: array

    field_vector: array
        subregion

    axisname: str
        'x', 'y', or 'z'

    outname: str
        output .h5 filename

    originalSize: float
        region of which we extracted this subset of fields from the original box, in code unit

    Returns
    -------

    """

    originalLevels = np.log2(1. / np.unique(dx_vector))
    numLevel = len(originalLevels)   # or numLevel = len(np.unique(dx_vector))
    levels = originalLevels               # let's just keep as
    originalLevels_vector = np.log2(1. / dx_vector)
    (originalLevels_vector == np.sort(originalLevels_vector)).all() == True

    levels_ii = originalLevels_vector
    levels_ii = np.array(levels_ii, dtype=int)
    levels = np.array(levels, dtype=int)

    d = {}
    for LL in levels:
        # find indices where they correspond to each of the levels
        mask = levels_ii == LL
        d[str(LL)] = mask
        print LL, np.max(levels_ii[mask]), np.min(levels_ii[mask])

    # check originalLevels_vector[ii] is properly saved to the right dict key
    ilowestLevel = d[str(int(levels.min()))]

    # coarset level
    Llow = str(int(levels.min()))
    pos = loc_vector[d[Llow], :]
    II = field_vector[d[Llow]]

    ilev = levels_ii[d[Llow]]
    assert ilev.min() == ilev.max() == int(Llow)

    print field_vector
    print II

    highestRes = 2.**(levels.max() * -1)
    N = originalSize / highestRes
    print N

    # make sure it's even number
    # N = 196
    N = int(N)
    if N % 2:  # odd number
        N += 1
    assert N % 2 == 0

    Ninit = N / len(levels) / 2
    field_cube = np.zeros((Ninit, Ninit, Ninit))

    imatrix = np.ones((2, 2, 2))

    for lll in np.sort(levels):
        if debug:
            print lll

        field_cube = np.kron(field_cube, imatrix)

        pos = loc_vector[d[str(lll)], :]

        xx = (pos[:, 0] - loc_vector[:, 0].min()) / \
            (loc_vector[:, 0].max() - loc_vector[:, 0].min())
        yy = (pos[:, 1] - loc_vector[:, 1].min()) / \
            (loc_vector[:, 1].max() - loc_vector[:, 1].min())
        zz = (pos[:, 2] - loc_vector[:, 2].min()) / \
            (loc_vector[:, 2].max() - loc_vector[:, 2].min())

        if debug:
            print '  x:', xx.min(), xx.max()
            print '  y:', yy.min(), yy.max()
            print '  z:', zz.min(), zz.max()
            #import pdb; pdb.set_trace()

        xx = xx * field_cube.shape[0]
        yy = yy * field_cube.shape[1]
        zz = zz * field_cube.shape[2]


        xpos = np.array(xx, dtype=int)
        ypos = np.array(yy, dtype=int)
        zpos = np.array(zz, dtype=int)

        xpos[xpos == field_cube.shape[0]] = field_cube.shape[0] - 1
        ypos[ypos == field_cube.shape[1]] = field_cube.shape[1] - 1
        zpos[zpos == field_cube.shape[2]] = field_cube.shape[2] - 1
        field_cube[xpos, ypos, zpos] = field_vector[d[str(lll)]]

        if debug:
            print '  field    :', field_cube.min(), field_cube.max()
            print '  field log:', np.log10(field_cube.min()), np.log10(field_cube.max())

    if debug:
        size = field_cube.shape[0]
        print "1d length of final resampled cube: ", size

        import matplotlib.pyplot as plt
        # 2**4: to chop off spurious edges
        plt.figure()
        plt.imshow(field_cube[:, :, 2**4:].sum(axis=0))
        plt.title('vel_' + axisname)
        plt.show()

    import h5py, os
    # if not exist, then create the .h5 file, else append mode
    if os.path.exists(outname):
        mode = "a"
    else:
        mode = "w"

    f = h5py.File(outname, mode)
    # create a dataset at the root
    f.create_dataset("/vel_" + axisname, data=field_cube)
    f.close()

    if debug:
        # read the h5 file back to make sure appended field
        f = h5py.File(outname, "r")
        print f.keys()
        f.close()

    return None



fieldsDict = {}
for i in fields:
    fieldsDict[i] = ds[i]

# print fieldsDict
# import pdb; pdb.set_trace()

for kkk, vvv in fieldsDict.iteritems():
    resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, debug=False)

# -- velocity --
vel = ds['vel']
axes = {'x': vel[:, 0], 'y': vel[:, 1], 'z': vel[:, 2]}

for kkk, vvv in axes.iteritems():
    resam_vel_field(dx_vector, loc_vector, vvv, kkk, outname, debug=True)


# ------