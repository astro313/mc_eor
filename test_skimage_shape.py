'''

script to test some skimage functions.. trying to understand what input shape do blob_dog and blob_log work on.. and whether we need first run rgb2gray()..

The resuts are confusing.. it seems like the functions can directly run on the cube, w/o needing to run thru rgb2gray() first, but results from which to trust..? They don't seem identical.

last mod: 6 July 2018


'''

print(__doc__)

from skimage.color import rgb2gray
from skimage.feature import blob_dog, blob_log
import matplotlib.pyplot as plt
import numpy as np
from skimage import data

blob2d = data.binary_blobs(length=128, n_dim=2)
blob3d = data.binary_blobs(length=128, n_dim=3)
print blob3d.shape
blob3d_gray = rgb2gray(blob3d)
print blob3d_gray.shape
# ok... after passing through rgb2gray, it collapsed into 2D....


# is running through rgb2gray() necessary before running skimage blob_dog and blob_log though?

# test dog 2d
blobs_dog = blob_dog(blob2d) # , max_sigma=10)
print blobs_dog
blobs_dog = blob_dog(rgb2gray(blob2d))   # , max_sigma=10)
print blobs_dog
# -- whether or not running rgb2gray() doesn't actually affect blob_dog() operation, at least on 2D image


# test dog 3d
blobs_dog = blob_dog(blob3d)
print blobs_dog
print blob3d.shape
blobs_dog = blob_dog(rgb2gray(blob3d))
print blobs_dog
print rgb2gray(blob3d).shape
# -- it does matter in the 3D case though... hmm..
# []
# (128, 128, 128)
# [[  40.        127.         16.777216]
#  [   0.         38.         16.777216]]
# (128, 128)



# repeat for LoG
blobs_log = blob_log(blob2d) # , max_sigma=10)
print blobs_log
print blobs_log.shape
blobs_log = blob_log(rgb2gray(blob2d))   # , max_sigma=10)
print blobs_log
print blobs_log.shape
# -- whether or not running rgb2gray() doesn't actually affect blob_log() operation, at least on 2D image


# test log 3d
blobs_log = blob_log(blob3d)
print blobs_log
print blobs_log.shape
print blob3d.shape
# (4170, 4)
# (128, 128, 128)
# ok, so it does work on the 3D cube.. and does output a (p,r,q,s) shape thing

blobs_log = blob_log(rgb2gray(blob3d))
print blobs_log
print blobs_log.shape
print rgb2gray(blob3d).shape
# -- so then, can blob_log work on the cube directly (without running through rgb2gray() first??)
# (99, 3)
# (128, 128)



# ---------
