'''

compare results from LoG and DoG


'''

from math import sqrt
from skimage import data
from skimage.feature import blob_dog, blob_log
from skimage.color import rgb2gray

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


_max = density.max()
image = density
image_gray = rgb2gray(image)

threshold = 100
factor = image_gray.max() / _max
threshold = 100.0 * factor
print threshold

blobs_log = blob_log(image_gray, max_sigma=30, num_sigma=10, threshold=threshold)

# Compute radii in the last column.
blobs_log[:, 3] = blobs_log[:, 3] * sqrt(3)


blobs_dog = blob_dog(image_gray, max_sigma=30, threshold=threshold)
blobs_dog[:, 3] = blobs_dog[:, 3] * sqrt(3)


blobs_list = [blobs_log, blobs_dog]
colors = ['yellow', 'lime']
titles = ['Laplacian of Gaussian', 'Difference of Gaussian']
sequence = zip(blobs_list, colors, titles)

fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)
ax = axes.ravel()

flatten_image = image.sum(axis=0)    # for imshow

for idx, (blobs, color, title) in enumerate(sequence):
    ax[idx].set_title(title)
    ax[idx].imshow(flatten_image, interpolation='nearest')
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
        ax[idx].add_patch(c)
    ax[idx].set_axis_off()

plt.tight_layout()
plt.show()