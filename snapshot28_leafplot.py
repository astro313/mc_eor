'''

plot physical properties of clouds.

last mod: 18 July 2018

'''


import cPickle as pickle
from snapshot_leafprop import Cloud
import numpy as np
import os
from plot_modules.plot_cloud_prop import setup_plot
from io_modules.manipulate_fetch_gal_fields import get_units, get_dx

plt = setup_plot()


# ----------- plot for a single snapshot, results from different ncuts ---

snapshot_num = 28
leafdir_out = "/mnt/home/daisyleung/mc_eor/test_brute/leaf_fields_" + str(snapshot_num) + '/'
if not os.path.isdir(leafdir_out):
    os.mkdir(leafdir_out)

from io_modules.leaf_pickle import load_in_pickled_leaf_singleSS
from plot_modules.plot_cloud_prop import unpack_xy, plot_stuff, plot_stuff_3dim

ss = load_in_pickled_leaf_singleSS(leafdir_out, snapshot_num)
to_plot = unpack_xy(ss)

plot_stuff("cloud mass", "mass over jeans mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff("cloud mass", "alpha vir", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff("gas sd", "sfr sd", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff("size pc", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff("gas sd", "sigma kms", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff("R2 pc2", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff_3dim("tff Myr", "size pc", "cloud mass", leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)
plot_stuff('gas sd cgs', 'sigmaSq over size', leglabel="ncut: ", to_plot=to_plot, outdir=leafdir_out)


from plot_modules.plot_cloud_prop import get_masses_all_clouds, massFuncUnbinnedCDF, massFuncPDF, massFuncDifferential
allmass = get_masses_all_clouds(ss)
tag = 'ss'+str(snapshot_num) +'diffncuts'
massFuncUnbinnedCDF(allmass, outdir=leafdir_out, tag=tag)
massFuncPDF(allmass, outdir=leafdir_out, tag=tag)
massFuncDifferential(allmass, outdir=leafdir_out, tag=tag)


# ---