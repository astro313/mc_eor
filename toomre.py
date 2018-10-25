'''

Toomre for gas component and stars and the effective combo of both

'''

import yt
import numpy as np
import cPickle as pickle
import pymses
from plot_modules.plot_cloud_prop import setup_plot
setup_plot()
from io_modules.manipulate_fetch_gal_fields import import_fetch_gal, prepare_unigrid, prepare_star_unigrid, get_units, import_fetch_stars

# load camera stuff
folder = 'precomputed_data/'
f_camera = folder + 'camera_settings.log'

with open(f_camera, 'rb') as f:
    cameraDat = pickle.load(f)

# prepare input data, if not already passed as arguments
isnap = 28
read_proper_unit = True

if read_proper_unit:
    import pymses
    ro = pymses.RamsesOutput("output", isnap)
    factor_vel = get_units(ro=ro)['vel'][0]
    factor_rho = get_units(ro=ro)['rho'][0]
    factor_R = get_units(ro=ro)['dx'][0] * 1.e-3
    region_size = cameraDat[str(isnap)]['size']
    region_size_kpc = region_size * factor_R
else:
    region_size_kpc = None

data = import_fetch_gal(isnap=isnap)
print data.keys()
ds, dd = prepare_unigrid(data=data,
                         add_unit=True,
                         regionsize_kpc=region_size_kpc,
                         debug=False)
# ds: yt StreamDataset
# dd: TYRegion
print dd['density'].shape
print dd['velx'].shape

# def _vel(field, data):
#     vel = np.c_[data['velx'], data['vely'], data['velz']]  # km/s
#     return vel

# ds.add_field(("velocity"), function=_vel)
# # _dd = ds.all_data()
# # print _dd['velocity']     # km/s

# project velocity by mass-weighting, project surface density
# axis : corresponding to the axis to slice along

axes = {'0': 'x', '1': 'y', '2': 'z'}
vel = {}
projected_totalGasSurfaceDensity = {}
coords = {}
center_plane = {}

xxx = dd["x"].reshape((224, 224, 224))
yyy = dd["y"].reshape((224, 224, 224))
zzz = dd["z"].reshape((224, 224, 224))
center = dd.get_field_parameter('center')

for kk, vv in axes.iteritems():

    proj = ds.proj('velx', int(kk), weight_field='h2density')
    velx_projected = proj['velx']
    proj = ds.proj('vely', int(kk), weight_field='h2density')
    vely_projected = proj['vely']
    proj = ds.proj('velz', int(kk), weight_field='h2density')
    velz_projected = proj['velz']
    # proj.save_as_dataset(filename='123_proj.h5')
    # proj_ds = yt.load("123_proj.h5")
    # p = yt.ProjectionPlot(proj_ds, int(kk), "velx")
    # p.save('123.pdf')

    if kk is '2':
        vel[kk] = np.c_[velx_projected, vely_projected]  # , velz_projected]
    elif kk is '0':
        vel[kk] = np.c_[vely_projected, velz_projected]  # , velx_projected]
    elif kk is '1':
        vel[kk] = np.c_[velx_projected, velz_projected]  # , vely_projected]

    # project total gas surface density
    proj = ds.proj("density", int(kk), method='integrate')
    projected_SD = proj['density']  # 1D-array
    # Does the unit here make sense? or do we need to multiply by some factor?

    projected_totalGasSurfaceDensity[kk] = projected_SD.reshape(
        (-1, int(np.sqrt(projected_SD.shape[0]))))
    # proj.save_as_dataset(filename='123_proj.h5')
    # proj_ds = yt.load("123_proj.h5")
    # p = yt.ProjectionPlot(proj_ds, "x", "density",
    #                       weight_field="density")
    # p.save('123.pdf')

    # project plane, coordinates.
    if kk is '2':
        x_slice = xxx[:, :, 0]
        y_slice = yyy[:, :, 0]
        coords[kk] = np.c_[x_slice.flatten(), y_slice.flatten()]
        center_plane[kk] = [center[0], center[1]]
    elif kk is '1':
        x_slice = xxx[:, 0, :]
        z_slice = zzz[:, 0, :]
        coords[kk] = np.c_[x_slice.flatten(), z_slice.flatten()]
        center_plane[kk] = [center[0], center[2]]
    elif kk is '0':
        y_slice = yyy[0, :, :]
        z_slice = zzz[0, :, :]
        coords[kk] = np.c_[y_slice.flatten(), z_slice.flatten()]
        center_plane[kk] = [center[1], center[2]]


plane = '2'
vel_plane = vel[plane]
projected_totalGasSurfaceDensity_plane = projected_totalGasSurfaceDensity[plane]
SD = projected_totalGasSurfaceDensity_plane.flatten()
print SD.max()
coords_plane = coords[plane]

# radial velocity
i_hat = coords_plane[:, 0] - center_plane[plane][0]
j_hat = coords_plane[:, 1] - center_plane[plane][1]
R = np.sqrt(i_hat**2 + j_hat**2)    # + k_hat**2)
i_hat /= R
j_hat /= R
# k_hat /= R
vi = vel_plane[:, 0]
vj = vel_plane[:, 1]
print vi
# vk = vel_plane[:, 2]
radial_vel = (vi * i_hat + vj * j_hat)   # + vk*k_hat)
radial_veloDisp = np.std(radial_vel)
print radial_vel
print radial_veloDisp

# v_phi
theta = np.arctan2(j_hat, i_hat)
_v_r = np.cos(theta) * vi + np.sin(theta) * vj
v_phi = (-np.sin(theta) * vi + np.cos(theta) * vi)
print v_phi
print np.std(v_phi)    # km/s

# kappa
# R is between 0 to 1
dv_phi_dr = np.gradient(v_phi, R)
kappa_sq = 2. * v_phi / R * (dv_phi_dr + v_phi / R)

# calculate Q_gas
A_gas = np.pi
A_stellar = 3.36
G = 6.67259e-8  # cgs

radial_veloDisp_cgs = radial_veloDisp * 1.e5
Q_gas = radial_veloDisp_cgs * np.sqrt(kappa_sq) / (A_gas * G * SD)
print Q_gas.shape

import matplotlib.pyplot as plt
plt.imshow(Q_gas.reshape((-1, 224)))
plt.savefig('123.pdf')

# stars
# read in particles
# read from resampled.h5
starData = import_fetch_stars(isnap=isnap, verbose=True, convert=True)   # default, convert=True
print starData.keys()

# max vel
# 2421836.4106407226 774885.1301048293 857806.5594254179 km/s, huhhh?! something wrong when regridding the velocities in cython_helper.pyx?!
ds, dd = prepare_star_unigrid(data=starData,
                              add_unit=True,     # since convert=True
                              regionsize_kpc=region_size_kpc,
                              debug=False)

# sigma, Kappa, Sigma


# Q_star =




# stars and gas effective Q, Romeo & Wiegert 2011
w = 2. * sigma_star * sigma_gas / (sigma_star**2 + sigma_gas**2)

WoverQstars = w / Q_star
WoverQgas = w / Q_gas

Q_twoComp_inv = np.empty((Q_star.shape[0], Q_star.shape[1]))

condition_1 = Q_star > Q_gas
condition_2 = Q_star < Q_gas

Q_twoComp_inv[condition_1] = WoverQstars + 1./Q_gas
Q_twoComp_inv[condition_2] = 1./Q_star + WoverQgas

Q_twoComp = 1./Q_twoComp_inv

print Q_twoComp.shape
plt.imshow(Q_twoComp.reshape((-1, 224)))
plt.savefig('789.pdf')



# Overplot clumps?