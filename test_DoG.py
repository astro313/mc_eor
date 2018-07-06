'''

Difference of Gaussian (DoG)

- faster approximation of LoG
- image is blurred with increasing standard deviations and the difference between two successively blurred images are stacked up in a cube.
- blobs are assumed to be light on dark background (white on black)

'''

import numpy as np
import matplotlib.pyplot as plt
from skimage.color import rgb2gray

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


from skimage.feature import blob_dog
# blobs are assumed to be light on dark background

# min_sigma = 1     # keep it low for smaller blobs
# max_sigma =       # keep it high to detect large blobs
# sigma_ratio =     # ratio btw the stddev of Gaussian Kernels used.
# overlap = 0.8     # between 0 and 1; if the area of two blobs overlaps by a fraction > this, smaller blob is eliminated.

# For each blob found, the method returns its coordinates and the standard deviation of the Gaussian kernel that detected the blob.
_max = density.max()
density = rgb2gray(density)
factor = density.max() / _max
threshold = 100.0 * factor
print threshold

csigma = blob_dog(density, threshold=threshold)
print csigma.shape
p = csigma[:, 0]
r = csigma[:, 1]
c = csigma[:, 2]
sigma = csigma[3]
radius_blob = np.sqrt(3) * sigma
# blob_log(density, min_sigma=min_sigma, max_sigma=max_sigma, )




#