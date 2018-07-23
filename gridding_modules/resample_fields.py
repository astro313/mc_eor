'''

resample AMR points to finest resolution.

Last mod: 13 July 2018


'''

import numpy as np

def resam_each_field(dx_vector, loc_vector, field_vector, fieldname, outname, originalSize, debug=True):

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

    levels_vec  = np.array(np.log2(1. / dx_vector),dtype = int)
    levels_uni  = np.unique(levels_vec)

    d = {}
    for LL in levels_uni:
        # find indices where they correspond to each of the levels
        mask = levels_vec == LL
        d[str(LL)] = mask
        print LL, np.max(levels_vec[mask]), np.min(levels_vec[mask])

    highestRes = 2.**(-levels_uni.max())
    N          = originalSize / highestRes
    if debug:
      print N

    # start with coarser level
    Ninit      = 1.*N / 2**len(levels_uni)
    Ninit      = int(np.ceil(Ninit))
    field_cube = np.zeros((Ninit, Ninit, Ninit))

    # set the kron magic multiplication
    imatrix    = np.ones((2, 2, 2))

    # actually this should be defined by the camera region
    bound_left  = np.min(loc_vector[:,:],axis = 0)
    bound_right = np.max(loc_vector[:,:],axis = 0)
    # below, the syntax as it appears in locate_parts
    '''
    # get the subregion
    center       = np.array(camera['center'])
    originalSize = np.array(camera['region_size'])
    if type(originalSize) == float is False:
        if len(originalSize) > 1:
            originalSize = originalSize[0]
    # init bound
    bound_left  = np.zeros((3))
    bound_right = np.zeros((3))
    # get the bounds
    bound_left   = center - originalSize/2.
    bound_right  = center + originalSize/2.
    '''

    if debug:
        aa = Ninit
        print 'init size',aa
        for lll in np.sort(levels_uni):
          aa = aa*2
        print 'final size',aa
        print levels_uni

    # cycle over levels
    for lll in np.sort(levels_uni):
        if debug:
            print 'level',lll

        # expand matrix with kron product
        field_cube = np.kron(field_cube, imatrix)

        # get cell position in the selected level
        pos        = loc_vector[d[str(lll)], :]
        # prepare cell id
        pos_id = np.zeros_like(pos)

        # cast positions to cube ids
        for iii in xrange(3):
          pos_id[:,iii] = (pos[:,iii] - bound_left[iii]) / (bound_right[iii] - bound_left[iii])
          pos_id[:,iii] = pos_id[:,iii] * field_cube.shape[iii]
        #
        pos_id = np.array(pos_id, dtype=int)
        # safety cut to avoid boundary problems
        for iii in xrange(3):
          pos_id[:,iii][pos_id[:,iii] >= field_cube.shape[iii]] = field_cube.shape[iii] - 1

        # assign field value to the cube
        field_cube[pos_id[:,0], pos_id[:,1], pos_id[:,2]] = field_vector[d[str(lll)]]

        if debug:
            print '  id        :'
            for iii in xrange(3):
              print '    ',iii,np.min(pos_id[:,iii]),np.max(pos_id[:,iii])
            print '  field     :', field_cube.min(), field_cube.max()
            print '  field log :', np.log10(field_cube.min()), np.log10(field_cube.max())


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
    import os
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

    return N

if __name__ == '__main__':

    # saved in fetch_gal_fields.py
    ssnum = 28
    ds = np.load('snapshot' + str(ssnum) + '_center_fields0123456-15.npz')
    outname = 'snapshot' + str(ssnum) + "_center_fields0123456-15_resampled.h5"

    import os
    import pickle

    folder = 'precomputed_data/'
    f_camera = folder + 'camera_settings.log'

    with open(f_camera,'rb') as f:
        data = pickle.load(f)

    originalSize = data[str(ssnum)]['size']

    import os
    if os.path.exists(outname):
        os.system('rm ' + outname)

    #
    dx_vector = ds['dx_vector']
    loc_vector = ds['loc_vector']

    fieldsToExtract = ['rho', 'P_nt', 'P', 'H2', 'Z']

    fieldsDict = {}
    for i in fieldsToExtract:
        fieldsDict[i] = ds[i]

    # print fieldsDict
    # import pdb; pdb.set_trace()

    for kkk, vvv in fieldsDict.iteritems():
      resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, originalSize, debug=True)

    # -- velocity --
    axes = {'vel_x': ds['vel'][:, 0], 'vel_y': ds['vel'][:, 1], 'vel_z': ds['vel'][:, 2]}

    for kkk, vvv in axes.iteritems():
      resam_each_field(dx_vector, loc_vector, vvv, kkk, outname, originalSize, debug=True)

# ------
