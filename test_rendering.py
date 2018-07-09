
import numpy as np
import matplotlib.pyplot as plt
import os
import yt
import h5py

Plot_stuff = False
debug      = False

fold_out   = 'test_png'

# convert from code unit density to g/cc (depending on how fetch_gal.py is implemented.)
convert_unit = True

if debug:
    namefile = "output/output_00028/info_00028.txt"
    myfile = namefile

    print "Loading file,", myfile
    pf = load(myfile)

f_in = 'snapshot28_center_fields012345-15_resampled.h5'

f       = h5py.File(f_in, "r")
density = f["rho"].value

if not os.path.isdir(fold_out):
    os.mkdir(fold_out)

if convert_unit:
    from fetch_gal_fields import get_units
    import pymses
    ro = pymses.RamsesOutput("output", 28)
    factor = get_units(ro = ro)['rho'][0]
    density *= factor
    print density.max()

data = dict(density=density)
ds = yt.load_uniform_grid(data, f["rho"].shape)
#dd = ds.all_data()

# from
# http://yt-project.org/doc/visualizing/volume_rendering.html
sc = yt.create_scene(ds, lens_type='perspective')
# Get a reference to the VolumeSource associated with this scene
# It is the first source associated with the scene, so we can refer to it
# using index 0.
source = sc[0]
# Set the bounds of the transfer function
source.tfh.set_bounds((1.e-2, 1.e+2))
# set that the transfer function should be evaluated in log space
source.tfh.set_log(True)

# Make underdense regions appear opaque
source.tfh.grey_opacity = True

# Plot the transfer function, along with the CDF of the density field to
# see how the transfer function corresponds to structure in the CDF
source.tfh.plot(fold_out+'/transfer_function.png', profile_field='density')

# save the image, flooring especially bright pixels for better contrast
sc.save(fold_out+'/rendering.png', sigma_clip=6.0)



