'''

2D contour finding...

Use a marching squares method to find constant valued contours in an image, array values are linearly interpolated to provide better precision of the output contours


http://scikit-image.org/docs/dev/auto_examples/edges/plot_contours.html#sphx-glr-auto-examples-edges-plot-contours-py



'''


import numpy as np
import matplotlib.pyplot as plt

from skimage import measure
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


def find2Dcontours(threeDdensity, axis=0, n_cut=100, step=100):
    '''

    Find 2D contours of constant levels using skimage package.

    Parameters
    ----------
    threeDdensity: 3D cube
        density cube

    axis: int
        along which axis to flatten in order for skimage to find 2D contours

    n_cut: float or int
        lowest level from which to start looking for contours

    step:
         multiplicative interval between contours

    N_cell_min: int
        not implemented at the moment
    '''

    axis = int(axis)

    c_min = n_cut
    c_max = 10**np.floor(np.log10(threeDdensity).max()+1)
    c_max = threeDdensity.max()

    print "min/max value for finding contours: ", c_min, c_max

    spacing = np.arange(c_min, c_max, step)

    # flatten image along x-axis
    flattened_image = threeDdensity.sum(axis=axis)


    # plot all contours found
    plt.figure(str(axis) + str(step))
    plt.imshow(np.log10(flattened_image), interpolation='nearest', cmap=plt.cm.gray)
    plt.colorbar()

    # Find contours at a constant value of spacing[ii]
    for ii in spacing:
        contours = measure.find_contours(flattened_image, ii)

        for n, contour in enumerate(contours):
            plt.plot(contour[:, 1], contour[:, 0], linewidth=2)

    plt.title('Collapse along ' + str(axis) + '-axis. Contours between steps: ' + str(step))
    plt.show(block=False)


find2Dcontours(density, 0, step=10)
find2Dcontours(density, 0, step=100)
find2Dcontours(density, 0, step=200)
find2Dcontours(density, 0, step=800)


find2Dcontours(density, 1, step=10)
find2Dcontours(density, 1, step=100)
find2Dcontours(density, 1, step=200)
find2Dcontours(density, 1, step=800)


find2Dcontours(density, 2, step=10)
find2Dcontours(density, 2, step=100)
find2Dcontours(density, 2, step=200)
find2Dcontours(density, 2, step=800)