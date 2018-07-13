'''

from http://yt-project.org/doc/visualizing/volume_rendering.html

Volume Render of H2 or H2 density or density, after resampled to finest grid.

Last mod: 10 July 2018


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

fold_out = 'test_png_28'

# convert from code unit density to g/cc (depending on how
# fetch_gal_fields.py is implemented.)
convert_unit = True

if debug:
    namefile = "output/output_00028/info_00028.txt"
    myfile = namefile

    print "Loading file,", myfile
    pf = load(myfile)

f_in = 'snapshot28_center_fields0123456-15_resampled.h5'

f          = h5py.File(f_in, "r")
density    = f["rho"].value
H2         = f["H2"].value
shape_data = f["rho"].shape

if not os.path.isdir(fold_out):
    os.mkdir(fold_out)

if convert_unit:
    from fetch_gal_fields import get_units
    import pymses
    ro = pymses.RamsesOutput("output", 28)
    factor = get_units(ro=ro)['rho'][0]      # 1/cm^3 (not H/cm^3)
    density *= factor
    print density.max()

from matplotlib import cm

if(1):
  # override for quick testing
  ris,__,__ = np.shape(H2)
  ris1,ris2 = ris/4 , 3*ris/4
  H2 = np.copy(density[ris1:ris2,ris1:ris2,ris1:ris2])
  H2[H2<=0] = np.min(H2[H2>0])
  shape_data = H2.shape

  bounds = (1.e-3, 1.e+3)

data = dict(density=H2)
# data = dict(density=density*H2)
#data = dict(density=density, H2=H2)
ds = yt.load_uniform_grid(data, shape_data)
#dd = ds.all_data()

# from refined_snapshot28_ytclump_fields import _h2density
# from yt.units import dimensions
# ds.add_field(("stream", "h2density"), function=_h2density, units='code_mass/code_length**3')   # units="1/cm**3")

sc = yt.create_scene(ds, lens_type='perspective')

# Get reference to the VolumeSource associated with this scene
# It is the first source associated with the scene, so we can refer to it
# using index 0.
source = sc[0]

source.set_log = True

# Set bounds of transfer function
# # plt.hist((density * H2).flatten(), bins=70)
# plt.hist((H2).flatten(), bins=70)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

if(0):
  # standard transfer function, see
  # http://yt-project.org/doc/visualizing/volume_rendering_tutorial.html#volume-rendering-tutorial
  tfh = TransferFunctionHelper(ds)
  tfh.set_field('density')
  tfh.set_log(True)
  tfh.set_bounds()
  tfh.build_transfer_function()
  tfh.tf.add_layers(5, colormap='inferno')

  # Grab the first render source and set it to use the new transfer function
  render_source = sc.get_source()
  render_source.transfer_function = tfh.tf
  #source.tfh.tf = tfh
  #source.tfh.bounds = bounds
elif(1):
  # select a continous distribution for values in the chosen range
  # transparency is linear
  # http://yt-project.org/doc/visualizing/volume_rendering.html
  tf = yt.ColorTransferFunction(np.log10(bounds))

  def linramp(vals, minval, maxval):
        return (vals - vals.min())/(vals.max() - vals.min())

  tf.map_to_colormap(np.log10(1.e-1), np.log10(1.e+3), colormap='inferno',scale_func=linramp)

  source.tfh.tf     = tf
  source.tfh.bounds = bounds

else:
  # transfer fucntion that select a discrete set of values in the chosen range
  # transparency is linear but it is awfu
  l
  # bounds = (1.E-01, 80.0)  # H2 density
  #bounds = (1.E-3, 0.1)      # H2 fraction
  tf = yt.ColorTransferFunction(np.log10(bounds))    #  Since this rendering is done in log space, the transfer function needs
  # to be specified in log space.

  # if H2 density:
  # tf.add_gaussian(np.log10(0.1), 0.001, [0.0, 0.0, 1.0, 1.0])
  # tf.add_gaussian(np.log10(10.0), 0.001, [0.0, 1.0, 0.0, 1.0])
  # tf.add_gaussian(np.log10(20.0), 0.001, [1.0, 0.0, 0.0, 1.0])

  # if H2 fraction:
  #tf.add_gaussian(np.log10(1.5e-3), 0.001, [0.0, 0.0, 1.0, 1.0])
  #tf.add_gaussian(np.log10(1.e-2), 0.001, [0.0, 1.0, 0.0, 1.0])
  #tf.add_gaussian(np.log10(0.07), 0.001, [1.0, 0.0, 0.0, 1.0])

  # density
  w         = 1.e-3
  #w         = 0.1
  cmap      = cm.get_cmap('viridis')   # meh
  list_cuts = [1.e-1,1.e+0,1.e+1,1.e+2,3.e+2]
  for i,v in enumerate(list_cuts):
    alp    = float(i+1)/len(list_cuts)
    print v
    print '  ',alp,np.log2(alp + 1),2**(alp - 1)
    #alp    = np.log2(alp + 1)
    #alp    = 2**(alp - 1)
    col    = list(cmap(alp))
    col[3] = alp
    print '  ',col
    tf.add_gaussian(np.log10(v), w , col)
  #tf.add_gaussian(np.log10(1.e-1), 0.001, [0.0, 0.0, 1.0, 0.2])
  #tf.add_gaussian(np.log10(1.e+0), 0.001, [0.0, 1.0, 0.0, 0.5])
  #tf.add_gaussian(np.log10(1.e+2), 0.001, [1.0, 0.0, 0.0, 1.0])
  source.tfh.tf     = tf
  source.tfh.bounds = bounds

  # Make underdense regions appear opaque
  source.tfh.grey_opacity = True

# Plot transfer function, along with the CDF of the density field to
# see how the transfer function corresponds to structure in the CDF
source.tfh.plot(fold_out + '/transfer_function.png', profile_field='density')

print sc.camera

if False:
  # change the camera, settings must be changed from pymses to yt internal definitions
  los_v = [0.85882976970482816, 0.49834986636750128, -0.11856996820546729]
  up_v  = [0.10255487134299716, 0.059509123032244614, 0.99294569974382518]
  cam   = sc.add_camera(ds, lens_type='perspective')
  cam.position = ds.arr([0.05,0.5,0.5],'code_length')
  cam.switch_orientation(normal_vector= los_v,
                             north_vector=up_v)

  cam.set_width(ds.domain_width * 1.0)
  print sc.camera
  #sc.render()

# save the image, flooring especially bright pixels for better contrast
#sc.save(fold_out + '/rendering.png', sigma_clip=6.0)
sc.save(fold_out + '/rendering.png', sigma_clip=4.0)
#sc.save(fold_out + '/rendering.png')

