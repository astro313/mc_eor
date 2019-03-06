'''

Toomre for star component.

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
cmap      = cm.get_cmap('viridis')

import astropy.constants as C
pc2cm = 3.086e18
cm2pc = 1/pc2cm
kg2Msun = 1. / 1.989E30
kg2g = 1.E+03
pc2cm = 3.086e18
J2erg   = 1.E+07
k_B_erg = C.k_B.value * J2erg
g2Msun = 1 / 1.989e33

# load camera stuff
folder = 'precomputed_data/'
f_camera = folder + 'camera_settings.log'

with open(f_camera, 'rb') as f:
    cameraDat = pickle.load(f)

# prepare input data, if not already passed as arguments
isnap = 28
read_proper_unit = True

verbose     = True
megaverbose = False

wg_var = 'mass'

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

# read in particles
# read from resampled.h5
starData = import_fetch_stars(
    isnap=isnap, verbose=verbose, convert=True)   # default, convert=True
if verbose:
  print starData.keys()

ds, dd = prepare_star_unigrid(data=starData,
                              add_unit=True,     # since convert=True
                              regionsize_kpc=region_size_kpc,
                              debug=False)

axes = {'0': 'x', '1': 'y', '2': 'z'}
vel = {}
veldisp = {}
veldisp_vertical = {}
projected_starSurfaceDensity = {}
coords = {}
coords3d = {}
center_plane = {}

n_bins = int(np.ceil(len(dd['mass'])**(1. / 3)))
xxx = dd["x"].convert_to_units('kpc').reshape((n_bins, n_bins, n_bins))
yyy = dd["y"].convert_to_units('kpc').reshape((n_bins, n_bins, n_bins))
zzz = dd["z"].convert_to_units('kpc').reshape((n_bins, n_bins, n_bins))
center = dd.get_field_parameter('center')

if verbose:
  print 'max/min mass',dd['mass'].max(),dd['mass'].min()
factor_mass, unit_mass = get_units(ro=ro)['mass']
# mass_code_unit = dd['mass'].value/factor_mass

# new field
dx = abs(xxx[0][1:] - xxx[1][:-1]).min().convert_to_units('pc')
dy = abs(yyy[0][1:] - yyy[1][:-1]).min().convert_to_units('pc')
dz = np.diff(zzz[0][0]).min().convert_to_units('pc')
rho_star = dd['mass'] / dx / dy / dz
if verbose:
  print 'rho_star in Msun/pc^3: '
  print '  max/min', rho_star.min(), rho_star.max()
if megaverbose:
  print rho_star

# print np.sum(rho_star, axis=-1, dtype='float64') * dx

# rho of star from mass
def _rho_star(field, data):
    kpc2cm = 1.e3 * 3.086e18
    rho = (data['mass'].convert_to_units('g')) / dx.convert_to_units('cm') / dy.convert_to_units('cm') / dz.convert_to_units('cm')
    return rho

if megaverbose:
  print dd['mass']
  print dx, dy, dz

ds.add_field(("rho_star"), function=_rho_star,
              units='g/cm**3')
_dd = ds.all_data()
if megaverbose:
  print _dd['rho_star']

for kk, vv in axes.iteritems():

    velocityx = dd['velx'].reshape((n_bins, n_bins, n_bins))
    velocityy = dd['vely'].reshape((n_bins, n_bins, n_bins))
    velocityz = dd['velz'].reshape((n_bins, n_bins, n_bins))

    # velocity weight by stellar mass intsead of H2 density here
    proj = ds.proj('velx', int(kk), weight_field=wg_var)
    velx_projected = proj['velx']
    velx_projected = velx_projected.reshape(
        (-1, int(np.sqrt(velx_projected.shape[0]))))
    proj = ds.proj('vely', int(kk), weight_field=wg_var)
    vely_projected = proj['vely']
    vely_projected = vely_projected.reshape(
        (-1, int(np.sqrt(vely_projected.shape[0]))))
    proj = ds.proj('velz', int(kk), weight_field=wg_var)
    velz_projected = proj['velz']
    velz_projected = velz_projected.reshape(
        (-1, int(np.sqrt(velz_projected.shape[0]))))

    # mass to surface density
    proj = ds.proj('rho_star', int(kk), method='integrate')
    projected_starSurfaceDensity[kk] = proj['rho_star'].reshape(
        (-1, int(np.sqrt(proj['rho_star'].shape[0]))))

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
        coords[kk] = [xxx[:, 0, :], zzz[:, 0, :]]
        veldisp[kk] = [velocityx, velocityz]
        veldisp_vertical[kk] = np.std(velocityy, axis=1)
        coords3d[kk] = [xxx, zzz]
        center_plane[kk] = [center[0], center[2]]
    elif kk is '0':
        vel[kk] = [vely_projected, velz_projected]
        veldisp[kk] = [velocityy, velocityz]
        veldisp_vertical[kk] = np.std(velocityx, axis=0)
        coords3d[kk] = [yyy, zzz]
        coords[kk] = [yyy[0, :, :], zzz[0, :, :]]
        center_plane[kk] = [center[1], center[2]]

plane = '2'

vel_plane = vel[plane]
veldisp_plane = veldisp[plane]
veldisp_vertical_plane = veldisp_vertical[plane]
coords_plane = coords[plane]
coords3d_plane = coords3d[plane]
projected_starSurfaceDensity_plane = projected_starSurfaceDensity[plane]
starSD = projected_starSurfaceDensity_plane

starSD = starSD.convert_to_units('Msun/pc**2')
if verbose:
  print 'max/min star SD in Msun/pc^2', np.max(starSD), np.min(starSD)


plt.figure()
im = plt.imshow(np.log10(veldisp_vertical_plane.value),
                origin='lower', cmap=cmap)
plt.title('vdisp vertical')
cbar = plt.colorbar(im)
cbar.set_label(r"$\sigma_z$ [km s$^{-1}$]", fontsize=16)
plt.show(block=False)


# radial velocity
i_hat = coords_plane[0] - center_plane[plane][0]
j_hat = coords_plane[1] - center_plane[plane][1]
R = np.sqrt(i_hat**2 + j_hat**2)    # + k_hat**2)
if verbose:
  print R.min()
i_hat /= R
j_hat /= R
#
vi = vel_plane[0]
vj = vel_plane[1]
radial_vel = (vi * i_hat + vj * j_hat)
radial_veloDisp = np.std(radial_vel)
if verbose:
  print 'radial velocity           ', np.max(radial_vel), np.min(radial_vel)
  print 'radial velocity dispersion', radial_veloDisp

# P_nt and P contribute to the velocity field in any case, as their divergence are sources in the Euler eq.
i_hat_3d = coords3d_plane[0] - center_plane[plane][0]
j_hat_3d = coords3d_plane[1] - center_plane[plane][1]
R_3d = np.sqrt(i_hat_3d**2 + j_hat_3d**2)
if megaverbose:
  print R_3d.min()
i_hat_3d /= R_3d
j_hat_3d /= R_3d

wg = dd[wg_var].reshape((n_bins, n_bins, n_bins))
v_r = veldisp_plane[0] * i_hat_3d + veldisp_plane[1] * j_hat_3d
mean_v_r = np.sum(wg * v_r, axis=int(plane))/np.sum(wg, axis=int(plane))

if plane == '0':
  tmp = mean_v_r[np.newaxis,:,:]
if plane == '1':
  tmp = mean_v_r[:,np.newaxis,:]
if plane == '2':
  tmp = mean_v_r[:,:,np.newaxis]

sigma_r = np.sqrt(np.sum(wg * (v_r - tmp)**2, axis=int(plane))/np.sum(wg, axis=int(plane)))

plt.figure()
plt.imshow(np.log10(sigma_r.value), cmap=cmap, origin='lower')
plt.title(r'$\sigma_r$ weighted by ' + wg_var)
plt.colorbar()
plt.show(block=False)
if verbose:
  print 'max/min radial vdisp in km/s (from velocity)', np.max(sigma_r), np.min(sigma_r)


# v_phi
theta = np.arctan2(j_hat, i_hat)
_v_r = np.cos(theta) * vi + np.sin(theta) * vj
v_phi = (-np.sin(theta) * vi + np.cos(theta) * vj)
if verbose:
  print 'maxmin v_phi    ', np.max(v_phi), np.min(v_phi)
  print 'std v_phi', np.std(v_phi)    # km/s

# kappa
R_cm = R.value * 1.e3 * pc2cm
R = R_cm

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
omega_mw_kms_kpc = 220. / 8
if verbose:
  print 'omega [km/s/kpc]'
  print '  mean',(omega_measured_standard_unit).mean()      # km/s/kpc
  print '  MW  ',omega_mw_kms_kpc

plt.figure()
plt.imshow(np.log10(omega_measured))
plt.title('log omega')
plt.colorbar()
plt.show(block=False)


nbins = 100
# Annular bin edges and centers
bins = np.linspace(0, 1, nbins) * r_slice.max()
bin_centers = bins[:-1] + (bins[1:] - bins[:-1]) / 2.

# Count how many pixels fall into each radial bin
hist, _ = np.histogram(r_slice, bins)
hist[hist == 0] = 1
if megaverbose:
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
plt.figure()
plt.imshow(domega_dr)
plt.colorbar()
plt.show(block=False)


rotation_frequency = omega_measured
kappa_sq = 2 * rotation_frequency / r_slice * (
    2 * r_slice * rotation_frequency + r_slice**2 * domega_dr)
kappa_sq[kappa_sq < 0] = np.min(kappa_sq[kappa_sq > 0])
kappa = np.sqrt(kappa_sq)
if verbose:
  print 'max min kappa',np.max(kappa),np.min(kappa)
if megaverbose:
  print kappa      # in the MW, kappa ~ omega, which is also true here

plt.figure()
plt.imshow(np.log10(kappa), cmap=cmap, origin='lower')
plt.title('kappa')
plt.colorbar()
plt.show(block=False)

A_star = 3.36
G = 6.67259e-8  # cgs

SD = starSD.convert_to_cgs()
if verbose:
  print 'starSD: '
  print '  max/min',np.max(SD),np.min(SD)
if megaverbose:
  print starSD
radial_veloDisp_cgs = sigma_r * 1.e5
whnzero = np.where(SD.value != 0)[0]
whnzero = whnzero.reshape((-1, int(np.sqrt(len(whnzero)))))
Q_star = np.zeros(SD.shape) * np.nan
Q_star[whnzero] = radial_veloDisp_cgs * kappa[whnzero] / \
    (A_star * G * SD[whnzero].value)

if verbose:
  print "Q_star: "
  print '  max/min',np.max(Q_star),np.min(Q_star)
if megaverbose:
  print Q_star
  print(np.isnan(Q_star) == True).any()

# define plot boundary kpc
_xmin =  (coords_plane[0] - center_plane[plane][0]).min()
_xmax =  (coords_plane[0] - center_plane[plane][0]).max()
_ymin =  (coords_plane[1] - center_plane[plane][1]).min()
_ymax =  (coords_plane[1] - center_plane[plane][1]).max()

plt.figure()
plt.imshow(np.log10(Q_star), origin='lower',
    extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
plt.title('Log Qstar without blanking')
cbar = plt.colorbar()
cbar.set_label(r"$\log{Q_{\rm star}}$")
plt.show(block=False)

plt.figure()
plt.imshow(np.log10(SD.value / cm2pc**2 * g2Msun), origin='lower', extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
plt.title(r'$\Sigma_{\rm star}$')
cbar = plt.colorbar()
cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]")
plt.show(block=False)
# seems two orders of magnitude too high?

plt.figure()
plt.imshow(np.log10(kappa * 3.086e+16), origin='lower', extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
plt.title(r'$\log{\kappa}$')
cbar = plt.colorbar()
cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]")
plt.show(block=False)


# show all quantities in one figure
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.10, right=0.90, hspace=0.1, wspace=0.25)
ax = plt.subplot(221)
im = ax.imshow(np.log10(SD.value / cm2pc**2 * g2Msun), origin='lower', \
                extent=(_xmin, _xmax, _ymin, _ymax), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]")

ax = plt.subplot(222)
im = ax.imshow(np.log10(sigma_r), origin='lower', extent=(_xmin, _xmax, _ymin, _ymax),
    cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\sigma}$ [km\,s$^{-1}$]")
plt.show(block=False)

ax = plt.subplot(223)
im = ax.imshow(np.log10(kappa * 3.086e+16), origin='lower', cmap=cmap,
                extent=(_xmin, _xmax, _ymin, _ymax))
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]")

ax = plt.subplot(224)
import matplotlib as mpl
cmap_div = cm.get_cmap('RdBu')         # divergent cmap
im = ax.imshow(np.log10(Q_star), origin='lower', \
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
cbar.set_label(r"$\log{Q_{\rm star}}$", fontsize=16)

# plt.tight_layout()
plt.show(block=False)


# zoomin and smooth
# only show central region of the plot (i.e., on the main galaxy)
xspacing = (_xmax - _xmin)/len(Q_star)
xruler = np.arange(_xmin, _xmax, xspacing)
rightBound = np.argmin(abs(xruler - 1.5))
leftBound = np.argmin(abs(xruler + 1.5))

yspacing = (_ymax - _ymin)/len(Q_star)
yruler = np.arange(_ymin, _ymax, yspacing)
topBound = np.argmin(abs(yruler - 1.5))
bottomBound = np.argmin(abs(yruler + 1.5))

fig = plt.figure(figsize=(9, 8))
fig.subplots_adjust(left=0.10, right=0.90, hspace=0.3, wspace=0.25)
ax = plt.subplot(221)
im = ax.imshow(gaussian_filter(np.log10(SD.value/ cm2pc**2 * g2Msun), sigma=0.8)[bottomBound: topBound, leftBound:rightBound], origin='lower', extent=(xruler[leftBound], xruler[rightBound], yruler[bottomBound], yruler[topBound]), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]")

ax = plt.subplot(222)
im = ax.imshow(gaussian_filter(np.log10(sigma_r.value), sigma=0.8)[bottomBound: topBound, leftBound:rightBound], origin='lower', extent=(xruler[leftBound], xruler[rightBound], yruler[bottomBound], yruler[topBound]), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\sigma}$ [km\,s$^{-1}$]")
plt.show(block=False)

ax = plt.subplot(223)
im = ax.imshow(gaussian_filter(np.log10(kappa * 3.086e+16), sigma=0.8)[bottomBound: topBound, leftBound:rightBound], origin='lower', extent=(xruler[leftBound], xruler[rightBound], yruler[bottomBound], yruler[topBound]), cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]", fontsize=16)
plt.xlabel('kpc', fontsize=16)
plt.ylabel('kpc', fontsize=16)

ax = plt.subplot(224)
im = ax.imshow(gaussian_filter(np.log10(Q_star), sigma=0.55)[bottomBound: topBound, leftBound:rightBound],
               origin='lower',
               extent=(xruler[leftBound],
                       xruler[rightBound],
                       yruler[bottomBound],
                       yruler[topBound]),
               cmap=cmap_div,
               vmin=-1, vmax=1     # clip at -1 < log10(Q_star) < 1
               )
cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                     ticks=[-1, 0, 1]
                    )
cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
cbar.set_label(r"$\log{Q_{\rm star}}$", fontsize=16)

# plt.tight_layout()
plt.show(block=False)

f_out = 'ss_' + str(isnap) + 'stars_toomre_proj_' + plane + '.png'
if verbose:
  print 'out to'
  print '  ',f_out
plt.savefig(f_out)


# #from scipy import interpolate
# f = interpolate.interp2d(xxx_gas, yyy_gas, Q_gas, kind='cubic')
# Q_gas_resampled = f(xxx, yyy)
# print Q_gas_resampled.shape
# assert Q_star.shape == Q_gas_resampled
# Q_gas = Q_gas_resampled

# # stars and gas effective Q, Romeo & Wiegert 2011
# w = 2. * radial_veloDisp_cgs * radial_veloDisp_gas_cgs.flatten() / (radial_veloDisp_cgs**2 + radial_veloDisp_gas_cgs.flatten()**2)

# WoverQstars = w / Q_star
# WoverQgas = w / Q_gas

# # 2D array
# Q_twoComp_inv = np.empty((Q_star.shape[0], Q_star.shape[1]))

# condition_1 = Q_star > Q_gas
# condition_2 = Q_star < Q_gas

# Q_twoComp_inv[condition_1] = WoverQstars + 1. / Q_gas
# Q_twoComp_inv[condition_2] = 1. / Q_star + WoverQgas

# Q_twoComp = 1. / Q_twoComp_inv
# print(np.isnan(Q_twoComp) == True).any()

# plt.figure()
# plt.imshow(np.log10(Q_twoComp).reshape((-1, Q_star.shape[0])))
# plt.show(block=False)





# def calculate_total_toomre_q(Q_star, Q_gas):
#     """Calculates the combined Toomre Q parameter

#     Uses a formula that takes into account the finite thickness of the disk
#     and that the gas and stars independently contribute to the gravitational
#     potential.

#     See Romeo & Wiegert (2011) [2011MNRAS.416.1191R] for details
#     http://adsabs.harvard.edu/abs/2011MNRAS.416.1191R

#     """
#     stars = self.stars
#     gas = self.gas

#     for container in [stars, gas]:
#         # This will implicitly find the velocity dispersions
#         if container.toomre_q is None:
#             self.calculate_toomre_q(container.field_type)

#     The effect of thickness is to increase the stability parameter of each
#     component by a factor T, which depends on the ratio of vertical to
#     radial velocity dispersion.
#     Eqn 8 of Romeo & Wiegert
#     T_s = 0.8 + 0.7 * (stars.velocity_dispersion_vertical /
#                        stars.velocity_dispersion_disk)
#     T_g = 0.8 + 0.7 * (gas.velocity_dispersion_vertical /
#                        gas.velocity_dispersion_disk)

#     sigma_g = gas.velocity_dispersion
#     sigma_s = stars.velocity_dispersion

#     W = 2 * sigma_s * sigma_g / (sigma_s**2 + sigma_g**2)

#     Q_s = Q_star
#     Q_g = Q_gas

      # Eqn 9 of Romeo & Wiegert
#     res1 = W / (Q_s * T_s) + 1 / (Q_g * T_g)
#     res2 = 1 / (Q_s * T_s) + W / (Q_g * T_g)

#     res = np.where(T_s * Q_s >= T_g * Q_g, res1, res2)

#     # Putting this on the gas component, totally ugly hack
#     self.gas.total_toomre_q = 1 / res

#     return res






