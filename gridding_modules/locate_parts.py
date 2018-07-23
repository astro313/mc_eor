'''

cast mass and age of star particles of finest resolution to unigrid cell resampled for AMR.

Last mod: 23 July 2018


'''

import numpy as np
import sys

sys.path.append('../')
from cy_modules.cython_helper import grid_particle_mass


def resample_star_mass_age(loc_vector, levels_vec, epoch_vector, mass_vector, \
                           outname, camera, Nrefined=None, debug=True):

    """

    cast each variable field of the star particles onto the same gridding as the AMR unigrid.

    Parameters
    ----------

    loc_vector: array

    levels_vec: array

    mass_vector: array
        mass of star particles in subregion

    epoch_vector: array
        epoch of star particles in subregion

    outname: str
        output .h5 filename

    camera: dict
        region of which we extracted this subset of fields from the original box, in code unit, with center and region_size

    Nrefined: int
       1D size of the resampled unigrid of the AMR stuff 

    debug: bool

    Returns
    -------

    """

    outfieldnames = ['mass', 'epoch']

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

    if Nrefined is None:
        dx_vector = np.array(1/(np.e * levels_vec))        
        levels_uni  = np.unique(levels_vec)

        highestRes = 2.**(-levels_uni.max())
        N = int(originalSize / highestRes)
    else:
        N = int(Nrefined)
    
    mass_cube = np.zeros((N, N, N))
    epoch_cube = np.zeros((N, N, N))

    xx = np.zeros_like(loc_vector)

    if debug:
      print 'max/min pos'
      for iii in xrange(3):
         print iii
         print '  pos  ',np.max(loc_vector[:,iii]),np.min(loc_vector[:,iii])
         print '  bound',bound_right[iii],bound_left[iii]

    for iii in xrange(3):
        xx[:,iii] = (loc_vector[:, iii] - bound_left[iii]) / \
                    abs(bound_right[iii] - bound_left[iii])

        xx[:,iii] = xx[:,iii] * mass_cube.shape[iii]

    xx        = np.round(xx,1)
    xx[xx>=N] = N-1
    xx[xx<=0]   = 0

    if debug:
      print 'max/min id'
      for iii in xrange(3):
         print iii,np.max(xx),np.min(xx)
    #xx = np.array(xx, dtype=int)

    mass_cube, epoch_cube = grid_particle_mass(mass_cube, epoch_cube, xx, mass_vector, epoch_vector)

    if debug:
        print '  Mass field    :', mass_cube.min(), mass_cube.max()
        print '  Mass field log:', np.log10(mass_cube[mass_cube>0].min()), np.log10(mass_cube.max())

        size = mass_cube.shape[0]
        print "1d length of final resampled cube: ", size


        to_plot = np.copy(mass_cube[:, :, :])
        print np.log10(to_plot[to_plot>0].min()), np.log10(to_plot.max())
        print to_plot.min(),to_plot.max()
        to_plot[to_plot<=0] = np.min(to_plot[to_plot>0])
        print to_plot.min(),to_plot.max()
        #import pdb; pdb.set_trace()

        import matplotlib.pyplot as plt
        # 2**4: to chop off spurious edges
        plt.imshow(np.log10(to_plot.sum(axis=0)))
        plt.title("Mass collapsed along one axis")
        plt.savefig("Mass_check.png")
        plt.show()

    import h5py
    import os
    # if not exist, then create the .h5 file, else append mode
    if os.path.exists(outname):
        mode = "a"
    else:
        mode = "w"

    print outname

    f = h5py.File(outname, mode)
    # create a dataset at the root
    f.create_dataset("/mass", data=mass_cube)
    f.create_dataset("/epoch", data=epoch_cube)
    f.close()

    if debug:
        # read the h5 file back to make sure appended field
        f = h5py.File(outname, "r")
        print f.keys()
        f.close()
    return None


if __name__ == '__main__':

    # saved in ../call_fetch_gal_fields.py
    ssnum = 28
    ds = np.load('../snapshot' + str(ssnum) + '_center_stars.npz')
    outname = '../snapshot' + str(ssnum) + "_center_stars_resampled.h5"

    import os
    import pickle

    if os.path.exists(outname):
        os.system('rm ' + outname)

    folder = '../precomputed_data/'
    f_camera = folder + 'camera_settings.log'

    with open(f_camera,'rb') as f:
        data = pickle.load(f)

    center      = data[str(ssnum)]['center_init']
    region_size = [data[str(ssnum)]['size']]          # subregion inside which we extracted, in code unit

    camera = {'center': center, 'region_size': region_size}

    loc_vec = ds['loc_vector']
    level_vec = ds['level']
    epoch_vec = ds['epoch']
    id_vec = ds['id']
    mass_vec = ds['mass']

    # hard-coded for testing
    amrfile = '../snapshot' + str(ssnum) + '_center_fields0123456-15_resampled.h5'
    import h5py
    f = h5py.File(amrfile, "r")
    N = f['rho'].shape[0]
                 
    # resample star particle mass and compute mass-weighted avg epoch correspondingly
    resample_star_mass_age(loc_vec, level_vec, epoch_vec, mass_vec, outname, camera, Nrefined=N, debug=True)

# ------
