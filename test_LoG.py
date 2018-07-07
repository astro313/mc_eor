'''

Laplacian of Gaussian (LoG)

- accurate but slow.
- computes the LoG images w/ successively increasing standard deviation and stacks them up in a cube. Blobs are local max in this cube.
- detecting larger blobs ususally slower because of larger kernel sizes during convolution.

'''


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
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


from skimage.feature import blob_log

# min_sigma = 1     # keep it low for smaller blobs
# max_sigma =       # keep it high to detect large blobs
# num_sigma = int(20)
# overlap = 0.8     # between 0 and 1; if the area of two blobs overlaps by a fraction > this, smaller blob is eliminated.

# For each blob found, the method returns its coordinates and the standard deviation of the Gaussian kernel that detected the blob.
_max = density.max()
# print 'shape before: ', density.shape
# density = rgb2gray(density)
# print 'shape after: ', density.shape
# factor = density.max() / _max
threshold = 100.0   # * factor
print threshold

print density.shape
csigma = blob_log(density, min_sigma=1, threshold=threshold)
print csigma.shape
p = csigma[:, 0]
r = csigma[:, 1]
c = csigma[:, 2]
sigma = csigma[:, 3]

print density.ndim
radius_blob = np.sqrt(3) * sigma
# blob_log(density, min_sigma=min_sigma, max_sigma=max_sigma, )

# visualize
flatten_image = density.sum(axis=0)    # for imshow
fig, ax = plt.subplots(1, 1)
ax.imshow(np.log10(flatten_image), interpolation='nearest')
for blob in csigma:
    y, x, c, r = blob
    cir = plt.Circle((c, x), r, color='blue', linewidth=2, fill=False)
    ax.add_patch(cir)
    ax.set_axis_off()
plt.tight_layout()
# plt.show()
plt.savefig('test_LoG_yz.png')


flatten_image = density.sum(axis=1)    # for imshow
fig, ax = plt.subplots(1, 1)
ax.imshow(np.log10(flatten_image), interpolation='nearest')
for blob in csigma:
    y, x, c, r = blob
    cir = plt.Circle((c, y), r, color='blue', linewidth=2, fill=False)
    ax.add_patch(cir)
    ax.set_axis_off()
plt.tight_layout()
# plt.show()
plt.savefig('test_LoG_xz.png')


flatten_image = density.sum(axis=2)    # for imshow
fig, ax = plt.subplots(1, 1)
ax.imshow(np.log10(flatten_image), interpolation='nearest')
for blob in csigma:
    y, x, c, r = blob
    cir = plt.Circle((x, y), r, color='blue', linewidth=2, fill=False)
    ax.add_patch(cir)
    ax.set_axis_off()
plt.tight_layout()
plt.show()
plt.savefig('test_LoG_xy.png')
#



#