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
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

from matplotlib import cm
cmap = cm.get_cmap('viridis')

import astropy.constants as C
kg2Msun = 1. / 1.989E30
kg2g = 1.E+03
pc2cm = 3.086e18
J2erg   = 1.E+07
k_B_erg = C.k_B.value * J2erg
g2Msun = 1 / 1.989e33
cm2pc = 1 / pc2cm


# load camera stuff
folder = 'precomputed_data/'
f_camera = folder + 'camera_settings.log'

with open(f_camera, 'rb') as f:
    cameraDat = pickle.load(f)

# prepare input data, if not already passed as arguments
isnap = 16
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

n_bins = int(np.ceil(len(dd['density'])**(1. / 3)))
xxx = dd["x"].reshape((n_bins, n_bins, n_bins))
yyy = dd["y"].reshape((n_bins, n_bins, n_bins))
zzz = dd["z"].reshape((n_bins, n_bins, n_bins))
center = dd.get_field_parameter('center')

# project velocity by mass-weighting, project surface density
# axis : corresponding to the axis to slice along

axes = {'0': 'x', '1': 'y', '2': 'z'}
vel = {}
veldisp = {}
veldisp_vertical = {}
projected_totalGasSurfaceDensity = {}
coords = {}
coords3d = {}
center_plane = {}


for kk, vv in axes.iteritems():

    velocityx = dd['velx'].reshape((n_bins, n_bins, n_bins))
    velocityy = dd['vely'].reshape((n_bins, n_bins, n_bins))
    velocityz = dd['velz'].reshape((n_bins, n_bins, n_bins))

    proj = ds.proj('velx', int(kk), weight_field='h2density')
    velx_projected = proj['velx']
    velx_projected = velx_projected.reshape(
        (-1, int(np.sqrt(velx_projected.shape[0]))))
    proj = ds.proj('vely', int(kk), weight_field='h2density')
    vely_projected = proj['vely']
    vely_projected = vely_projected.reshape(
        (-1, int(np.sqrt(vely_projected.shape[0]))))
    proj = ds.proj('velz', int(kk), weight_field='h2density')
    velz_projected = proj['velz']
    velz_projected = velz_projected.reshape(
        (-1, int(np.sqrt(velz_projected.shape[0]))))

    # project total gas surface density
    proj = ds.proj("density", int(kk), method='integrate')
    projected_totalGasSurfaceDensity[kk] = proj['density'].reshape(
        (-1, int(np.sqrt(proj['density'].shape[0]))))

    # project plane, coordinates.
    if kk is '2':
        vel[kk] = [velx_projected, vely_projected]
        veldisp[kk] = [velocityx, velocityy]
        veldisp_vertical[kk] = np.std(velocityz, axis=2)
        coords[kk] = [xxx[:, :, 0], yyy[:, :, 0]]
        coords3d[kk] = [xxx, yyy]
        center_plane[kk] = [center[0], center[1]]
    elif kk is '1':
        vel[kk] = [velx_projected, velz_projected]
        veldisp[kk] = [velocityx, velocityz]
        veldisp_vertical[kk] = np.std(velocityy, axis=1)
        coords[kk] = [xxx[:, 0, :], zzz[:, 0, :]]
        coords3d[kk] = [xxx, zzz]
        center_plane[kk] = [center[0], center[2]]
    elif kk is '0':
        vel[kk] = [vely_projected, velz_projected]
        veldisp[kk] = [velocityy, velocityz]
        veldisp_vertical[kk] = np.std(velocityx, axis=0)
        coords3d[kk] = [yyy, zzz]
        coords[kk] = [yyy[0, :, :], zzz[0, :, :]]
        center_plane[kk] = [center[1], center[2]]


plane = '0'

vel_plane = vel[plane]
veldisp_plane = veldisp[plane]
veldisp_vertical_plane = veldisp_vertical[plane]
coords_plane = coords[plane]
coords3d_plane = coords3d[plane]
projected_totalGasSurfaceDensity_plane = projected_totalGasSurfaceDensity[plane]
SD = projected_totalGasSurfaceDensity_plane
print 'max/min SD', np.max(SD), np.min(SD)


plt.figure()
im = plt.imshow(veldisp_vertical_plane.value, origin='lower', cmap=cmap)
plt.title('vdisp vertical')
cbar = plt.colorbar(im)
cbar.set_label(r"$\sigma_z$ [km s$^{-1}$]")
plt.show(block=False)


# radial velocity
i_hat = coords_plane[0] - center_plane[plane][0]
j_hat = coords_plane[1] - center_plane[plane][1]
R = np.sqrt(i_hat**2 + j_hat**2)    # + k_hat**2)
print R.min()
i_hat /= R
j_hat /= R
#
vi = vel_plane[0]
vj = vel_plane[1]
radial_vel = (vi * i_hat + vj * j_hat)
radial_veloDisp = np.std(radial_vel)
print 'radial velocity           ', np.max(radial_vel), np.min(radial_vel)
print 'radial velocity dispersion', radial_veloDisp

# P_nt and P contribute to the velocity field in any case, as their divergence are sources in the Euler eq.
i_hat_3d = coords3d_plane[0] - center_plane[plane][0]
j_hat_3d = coords3d_plane[1] - center_plane[plane][1]
R_3d = np.sqrt(i_hat_3d**2 + j_hat_3d**2)
print R_3d.min()
i_hat_3d /= R_3d
j_hat_3d /= R_3d

# weigh by mass
wg = dd['h2density'].reshape((n_bins, n_bins, n_bins))
v_r = veldisp_plane[0] * i_hat_3d + veldisp_plane[1] * j_hat_3d
mean_v_r = np.sum(wg * v_r, axis=2)/np.sum(wg, axis=2)
sigma_r = np.sqrt(np.sum(wg * (v_r - mean_v_r[:,:,np.newaxis])**2, axis=2)/np.sum(wg, axis=2))
plt.figure()
plt.imshow(sigma_r.value, cmap=cmap, origin='lower')
plt.title(r'$\sigma_r$ weighted by h2 density')
plt.colorbar()
plt.show(block=False)
print 'max/min radial vdisp in km/s (from velocity)', np.max(sigma_r), np.min(sigma_r)


_sigma_r = np.std(veldisp_plane[0] * i_hat + veldisp_plane[1] * j_hat, axis=2)
plt.figure()
plt.imshow(_sigma_r.value, cmap=cmap, origin='lower')
plt.title(r'$\sigma_r$ not weighted')
plt.colorbar()
plt.show(block=False)
print 'max/min radial vdisp in km/s (from velocity)', np.max(_sigma_r), np.min(_sigma_r)


# v_phi
theta = np.arctan2(j_hat, i_hat)
_v_r = np.cos(theta) * vi + np.sin(theta) * vj
v_phi = (-np.sin(theta) * vi + np.cos(theta) * vj)
print 'maxmin v_phi    ', np.max(v_phi), np.min(v_phi)
print 'std v_phi', np.std(v_phi)    # km/s

# kappa
R_cm = R.value * 1.e3 * pc2cm
R = R_cm

# smooth map to regularized the derivates
sigma_r = gaussian_filter(sigma_r, 0.15)
v_phi = gaussian_filter(v_phi, 0.15)
SD = gaussian_filter(SD, 0.15)

try:
  v_phi = v_phi.value * 1.e5         # cm/s
except AttributeError:
  v_phi = v_phi * 1.e5

cost = i_hat
sint = j_hat
x_slice = (coords_plane[0] - center_plane[plane][0]) * 1.e3 * pc2cm
x_slice = x_slice.value
y_slice = (coords_plane[1] - center_plane[plane][1]) * 1.e3 * pc2cm
y_slice = y_slice.value
r_slice = np.sqrt(x_slice**2 + y_slice**2)    # cm


omega_measured = v_phi / r_slice
omega_measured_standard_unit = v_phi / 1.e5 / (r_slice / pc2cm / 1.e3)
print(omega_measured_standard_unit).mean()      # km/s/kpc
omega_mw_kms_kpc = 220. / 8
print omega_mw_kms_kpc
# ok, unit comparable to MW value

# plt.figure()
# plt.imshow(omega_measured)
# plt.colorbar()
# plt.show(block=False)


nbins = 100
# Annular bin edges and centers
bins = np.linspace(0, 1, nbins) * r_slice.max()
bin_centers = bins[:-1] + (bins[1:] - bins[:-1]) / 2.

# Count how many pixels fall into each radial bin
hist, _ = np.histogram(r_slice, bins)
hist[hist == 0] = 1
print hist

# Get flat 1D indices into r_slice for each bin.
inds = np.digitize(r_slice.flat, bins) - 1

# Calculate mean omega at each radius.
# Need to append [:-1] to get rid of counts outside bin range.
omega_of_r = np.bincount(inds, weights=omega_measured.flat)[:-1] / hist

# Calculate standard deviation of the mean omega for each bin.
omega2_of_r = np.bincount(inds, weights=omega_measured.flat[:]**2)[:-1] / hist
omega_of_r_std = np.sqrt(omega2_of_r - omega_of_r**2)

# Calculate radial derivative of the rotation frequency using a spline
# interpolator
omega = np.zeros(r_slice.shape)
omega_deriv = np.zeros(r_slice.shape)

from scipy import interpolate
omega_interpolator = interpolate.splrep(
    bin_centers, omega_of_r, k=5, w=1 / omega_of_r_std)

omega.flat[:] = interpolate.splev(
    r_slice.flat, omega_interpolator)
rotation_frequency = omega       # should be 1/s

omega_deriv.flat[:] = interpolate.splev(
    r_slice.flat, omega_interpolator, der=1)
rotation_frequency_derivative = omega_deriv

domega_dr = np.zeros(r_slice.shape)
domega_dr.flat[:] = interpolate.splev(
    r_slice.flat, omega_interpolator, der=1)

# Finally, calculate epicyclic frequency.
# plt.figure()
# plt.imshow(domega_dr)
# plt.colorbar()
# plt.show(block=False)


rotation_frequency = omega_measured
kappa_sq = 2 * rotation_frequency / r_slice * (
    2 * r_slice * rotation_frequency + r_slice**2 * domega_dr)
kappa_sq[kappa_sq < 0] = np.min(kappa_sq[kappa_sq > 0])
kappa = np.sqrt(kappa_sq)
print kappa      # in the MW, kappa ~ omega, which is also true here

plt.figure()
plt.imshow(np.log10(kappa), cmap=cmap, origin='lower')
plt.colorbar()
plt.show(block=False)

# calculate Q_gas
A_gas = np.pi
A_stellar = 3.36
G = 6.67259e-8  # cgs

radial_veloDisp_cgs = sigma_r * 1.e5
whnzero = np.where(SD != 0)
Q_gas = np.zeros(SD.shape) * np.nan
Q_gas[whnzero] = radial_veloDisp_cgs[whnzero] * kappa[whnzero] / \
    (A_gas * G * SD[whnzero])

print "Q_gas: ", Q_gas
print(np.isnan(Q_gas) == True).any()
# qq = Q_gas[~np.isnan(Q_gas)]
# print qq[qq > 0].min()

# define plot boundary kpc
_xmin =  (coords_plane[0] - center_plane[plane][0]).min()
_xmax =  (coords_plane[0] - center_plane[plane][0]).max()
_ymin =  (coords_plane[1] - center_plane[plane][1]).min()
_ymax =  (coords_plane[1] - center_plane[plane][1]).max()

plt.figure()
plt.imshow(np.log10(Q_gas), origin='lower',
            extent=(_xmin, _xmax, _ymin, _ymax)
            #, cmap=cmap
            )
plt.title('Log Qgas without blanking')
cbar = plt.colorbar()
cbar.set_label(r"$\log{Q_{\rm gas}}$")
plt.show(block=False)

# show all quantities in one figure
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.10, right=0.90, hspace=0.1, wspace=0.25)
ax = plt.subplot(221)
im = ax.imshow(np.log10(SD / cm2pc**2 * g2Msun), origin='lower', extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]")

ax = plt.subplot(222)
im = ax.imshow(sigma_r, origin='lower', extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\sigma$ [km\,s$^{-1}$]")
plt.show(block=False)

ax = plt.subplot(223)
im = ax.imshow(np.log10(kappa * 3.086e+16), origin='lower', extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]")

ax = plt.subplot(224)
import matplotlib as mpl
cmap_div = cm.get_cmap('RdBu')         # divergent cmap
im = ax.imshow(np.log10(Q_gas), origin='lower', \
               extent=(_xmin, _xmax, _ymin, _ymax),
               cmap=cmap_div,
               vmin=-1,
               vmax=1
               )
cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                     ticks=[-1, 0, 1]
                    )
# cbar.set_clim(-1, 1)                      # set color max min range
cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
cbar.set_label(r"$\log{Q_{\rm gas}}$")

# plt.tight_layout()
plt.show(block=False)


# zoomin
# only show central region of the plot (i.e., on the main galaxy)
# from -1.5 kpc to 1.5 kpc
xspacing = (_xmax - _xmin)/len(Q_gas)
xruler = np.arange(_xmin, _xmax, xspacing)
rightBound = np.argmin(abs(xruler - 1.5))
leftBound = np.argmin(abs(xruler + 1.5))

yspacing = (_ymax - _ymin)/len(Q_gas)
yruler = np.arange(_ymin, _ymax, yspacing)
topBound = np.argmin(abs(yruler - 1.5))
bottomBound = np.argmin(abs(yruler + 1.5))

fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.10, right=0.90, hspace=0.25, wspace=0.25)
ax = plt.subplot(221)
im = ax.imshow(np.log10(SD / cm2pc**2 * g2Msun)[bottomBound: topBound, leftBound:rightBound],
              origin='lower',
              extent=(xruler[leftBound],
                      xruler[rightBound],
                      yruler[bottomBound],
                      yruler[topBound]),
              cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]")

ax = plt.subplot(222)
im = ax.imshow(sigma_r[bottomBound: topBound, leftBound:rightBound],
               origin='lower',
               extent=(xruler[leftBound],
                       xruler[rightBound],
                       yruler[bottomBound],
                       yruler[topBound]),
               cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\sigma$ [km\,s$^{-1}$]")
plt.show(block=False)

ax = plt.subplot(223)
im = ax.imshow(np.log10(kappa * 3.086e+16)[bottomBound: topBound, leftBound:rightBound],
               origin='lower',
               extent=(xruler[leftBound],
                       xruler[rightBound],
                       yruler[bottomBound],
                       yruler[topBound]),
               cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]")
plt.xlabel('kpc')
plt.ylabel('kpc')

ax = plt.subplot(224)
im = ax.imshow(np.log10(Q_gas)[bottomBound: topBound, leftBound:rightBound],
               origin='lower',
               extent=(xruler[leftBound],
                       xruler[rightBound],
                       yruler[bottomBound],
                       yruler[topBound]),
               cmap=cmap,
               vmin=-1, vmax=1     # clip at -1 < log10(Q_gas) < 1
               )
cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                     ticks=[-1, 0, 1]
                    )
cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
cbar.set_label(r"$\log{Q_{\rm gas}}$")

# plt.tight_layout()
plt.show(block=False)
# plt.savefig('ss_' + str(isnap) + '_toomre_proj_' + plane + '.png')


