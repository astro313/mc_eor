'''

Determinant of Hessian (DoH) - only for testing. Clumps identified are pretty coarse and not meaningful for us..

- fast, but only works on 2D images.
- detects blobs by finding maximas in the matrix of the Determinant of Hessian of the image
-  detection speed is independent of the size of blobs as internally the implementation uses box filters instead of convolutions
- Bright on dark as well as dark on bright blobs are detected
- small blobs (<3px) are NOT detected accurately


Last mod: 5 July 2018



'''

import numpy as np
import matplotlib.pyplot as plt

import h5py

f = h5py.File("snapshot28_center_densityfield_resampled.h5", "r")
density = f["density"].value

# convert from code unit density to g/cc (depending on how fetch_gal.py is implemented.)
convert_unit = True

if convert_unit:
    import pymses
    from pymses.utils import constants as C

    ro = pymses.RamsesOutput("output", 28)
    factor = ro.info["unit_density"].express(C.H_cc)
    density *= factor
    print density.max()


def find2Dcontours(threeDdensity, axis=0, n_cut=100):
    '''

    blobs are assumed to be light on dark background

    each blob found, the method returns its coordinates and the standard deviation of the Gaussian kernel that detected the blob.


    Parameters
    ----------
    threeDdensity: 3D cube
        density cube

    axis: int
        along which axis to flatten in order for skimage to find 2D contours

    n_cut: float or int
        lowest level from which to start looking for contours

    N_cell_min: int
        not implemented at the moment

    overlap: float
        not implemented; between 0 and 1; if the area of two blobs overlaps by a fraction > this, smaller blob is eliminated.

    # min_sigma = 1     # keep it low for smaller blobs
    # max_sigma =       # keep it high to detect large blobs
    # num_sigma =

    '''

    axis = int(axis)

    from skimage.feature import blob_doh
    from skimage.color import rgb2gray

    c_min = n_cut

    # flatten image along some axis
    flattened_image = threeDdensity.sum(axis=axis)
    flattened_image = rgb2gray(flattened_image)

    # plot all contours found
    fig, ax = plt.subplots()
    ax.imshow(np.log10(flattened_image),
               interpolation='nearest', cmap=plt.cm.gray)
    coord_sigma = blob_doh(flattened_image, threshold=n_cut)

    for blob in coord_sigma:
        y, x, r = blob
        c = plt.Circle((x, y), r, color='r', linewidth=2, fill=False)
        ax.add_patch(c)
    plt.show(block=False)

    # blob_log(density, min_sigma=min_sigma, max_sigma=max_sigma, )

    return coord_sigma, flattened_image


csigma, im = find2Dcontours(density, axis=0)
y = csigma[0]
x = csigma[1]
sigma = csigma[2]
radius_blob = sigma


csigma, im = find2Dcontours(density, axis=1)
y = csigma[0]
x = csigma[1]
sigma = csigma[2]
radius_blob = sigma


csigma, im = find2Dcontours(density, axis=2)
y = csigma[0]
x = csigma[1]
sigma = csigma[2]
radius_blob = sigma


#