
import os
import pickle

this_folder = "./"
f_camera    = this_folder+'camera_settings.log'

with open(f_camera,'rb') as f:
  data = pickle.load(f)

for nout in sorted(data.keys(), key = lambda iout: int(iout)):
  cam = data[nout]
  print nout
  for key in sorted(cam.keys()):
    print '  ',key,cam[key]

