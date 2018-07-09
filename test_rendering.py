'''

from http://yt-project.org/doc/visualizing/volume_rendering.html

Volume Render of H2 or H2 density or density, after resampled to finest grid.

Last mod: 9 July 2018


'''

print(__doc__)

import numpy as np
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
import os
import yt
import h5py
import matplotlib.pyplot as plt

Plot_stuff = False
debug = False

fold_out = 'test_png'

# convert from code unit density to g/cc (depending on how
# fetch_gal_fields.py is implemented.)
convert_unit = True

if debug:
    namefile = "output/output_00028/info_00028.txt"
    myfile = namefile

    print "Loading file,", myfile
    pf = load(myfile)

f_in = 'snapshot28_center_fields012345-15_resampled.h5'

f = h5py.File(f_in, "r")
density = f["rho"].value
H2 = f["H2"].value

if not os.path.isdir(fold_out):
    os.mkdir(fold_out)

if convert_unit:
    from fetch_gal_fields import get_units
    import pymses
    ro = pymses.RamsesOutput("output", 28)
    factor = get_units(ro=ro)['rho'][0]      # 1/cm^3 (not H/cm^3)
    density *= factor
    print density.max()

data = dict(density=H2)
# data = dict(density=density*H2)
#data = dict(density=density, H2=H2)
ds = yt.load_uniform_grid(data, f["rho"].shape)
dd = ds.all_data()

# from refined_snapshot28_ytclump_fields import _h2density
# from yt.units import dimensions
# ds.add_field(("stream", "h2density"), function=_h2density, units='code_mass/code_length**3')   # units="1/cm**3")

sc = yt.create_scene(ds, lens_type='perspective')

# Get reference to the VolumeSource associated with this scene
# It is the first source associated with the scene, so we can refer to it
# using index 0.
source = sc[0]

# Set bounds of transfer function
# plt.hist((density * H2).flatten(), bins=70)
# plt.hist((H2).flatten(), bins=70)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# bounds = (1.E-01, 80.0)  # H2 density
bounds = (1.E-3, 0.1)
tf = yt.ColorTransferFunction(np.log10(bounds))    #  Since this rendering is done in log space, the transfer function needs
# to be specified in log space.

# H2 density
# tf.add_gaussian(np.log10(0.1), 0.001, [0.0, 0.0, 1.0, 1.0])
# tf.add_gaussian(np.log10(10.0), 0.001, [0.0, 1.0, 0.0, 1.0])
# tf.add_gaussian(np.log10(20.0), 0.001, [1.0, 0.0, 0.0, 1.0])

# H2 fraction
tf.add_gaussian(np.log10(1.5e-3), 0.001, [0.0, 0.0, 1.0, 1.0])
tf.add_gaussian(np.log10(1.e-2), 0.001, [0.0, 1.0, 0.0, 1.0])
tf.add_gaussian(np.log10(0.07), 0.001, [1.0, 0.0, 0.0, 1.0])

source.tfh.tf = tf
source.tfh.bounds = bounds

# Make underdense regions appear opaque
source.tfh.grey_opacity = True


# Plot transfer function, along with the CDF of the density field to
# see how the transfer function corresponds to structure in the CDF
source.tfh.plot(fold_out + '/transfer_function.png', profile_field='density')

# save the image, flooring especially bright pixels for better contrast
sc.save(fold_out + '/rendering.png', sigma_clip=4.0)

