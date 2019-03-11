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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter

import os,sys

from matplotlib import cm
cmap     = cm.get_cmap('viridis')
cmap2    = cm.get_cmap('magma')
cmap_div = cm.get_cmap('RdBu')         # divergent cmap

from pymses.utils import constants as C_py
from yt import units as C_yt
import astropy.constants as C_ap
import astropy.units as U


class ToomreAnalyze(object):
  """
  Single Toomre Q object
  equation taken from Inoue+16
     http://adsabs.harvard.edu/abs/2016MNRAS.456.2052I

  """
  def __init__(self, isnap, wg_var, field_type, plane, smooth_size_kpc=0.1, read_proper_unit=True, verbose=True, debug=False, convertPart=True, megaverbose = False, min_wg = 'min', fold_out = '' , show = False):

    self.isnap = isnap
    self.read_proper_unit = read_proper_unit
    #
    self.wg_var = wg_var         # density for gas, mass for stars
    self.min_wg = min_wg         # clipping method for the weight
    #
    self.debug = debug
    self.verbose = verbose
    self.show    = show
    self.megaverbose = megaverbose
    self.convertPart = convertPart
    self.plane = plane
    self.field_type = field_type
    #
    if self.field_type == 'gas':
      self.c_clump  = 'darkgoldenrod'
    elif self.field_type == 'star':
      self.c_clump = '#bcbd22'    # https://matplotlib.org/users/dflt_style_changes.html
    #
    self.smooth_size = 0.                  # size for smoothing image in code units
    self.smooth_kpc  = smooth_size_kpc     # size for smoothing image in kpc
    self.smooth_log  = False # i do think we should smooth in linear space

    assert self.plane in ['0', '1', '2']
    assert self.field_type in ['star', 'gas']

    self.kg2Msun  = 1. / 1.989E30
    self.kg2g     = 1.E+03
    self.pc2cm    = 3.086e18
    self.J2erg    = 1.E+07
    self.k_B_erg  = C_ap.k_B.value * self.J2erg
    self.g2Msun   = 1 / 1.989e33
    self.cm2pc    = 1 / self.pc2cm
    self.kpc2cm = 1.e3 * self.pc2cm

    self.A_gas = np.pi
    self.A_star = 3.36
    self.G = 6.67259e-8  # cgs

    self.cameraFolder = 'precomputed_data/'
    self.f_camera = self.cameraFolder + 'camera_settings.log'

    # project velocity by mass-weighting, project surface density
    # axis : corresponding to the axis to slice along'none'

    self.axes = {'0': 'x', '1': 'y', '2': 'z'}
    self.vel = {}
    self.veldisp = {}
    self.veldisp_vertical = {}
    self.c_s_eff_proj = {}
    self.projected_SurfaceDensity = {}
    self.coords = {}
    self.coords3d = {}
    self.center_plane = {}

    # regulate output
    self.fold_out       = fold_out
    if (self.fold_out !='' and not os.path.isdir(self.fold_out)):
      os.mkdir(self.fold_out)

  def load_cam_stuff(self):
    with open(self.f_camera, 'rb') as f:
        self.cameraDat = pickle.load(f)

    if self.read_proper_unit:
        import pymses
        self.ro = pymses.RamsesOutput("output", self.isnap)
        self.factor_vel = get_units(ro=self.ro)['vel'][0]
        self.factor_rho = get_units(ro=self.ro)['rho'][0]
        self.factor_R = get_units(ro=self.ro)['dx'][0] * 1.e-3
        self.region_size = self.cameraDat[str(self.isnap)]['size']
        self.region_size_kpc = self.region_size * self.factor_R
    else:
        self.region_size_kpc = None


  def load_data(self):

    if self.field_type == 'gas':
      self.data = import_fetch_gal(isnap=self.isnap
                                  ,clipping = self.min_wg
                                  )
      if self.megaverbose:
        print self.data.keys()
      self.ds, self.dd = prepare_unigrid(data=self.data,
                               add_unit=True,
                               regionsize_kpc=self.region_size_kpc,
                               debug=self.debug
                               )

    elif self.field_type == 'star':
      # read from resampled.h5
      self.starData = import_fetch_stars(isnap=self.isnap,
                                         verbose=self.verbose,
                                         convert=self.convertPart, clipping = self.min_wg)
      if self.megaverbose:
        print self.starData.keys()
      self.ds, self.dd = prepare_star_unigrid(data=self.starData,
                                    add_unit=True,     # since convert=True
                                    regionsize_kpc=self.region_size_kpc,
                                    debug=self.debug)
    return self.ds, self.dd


  def setup_config(self):
    self.n_bins = int(np.ceil(len(self.dd[self.wg_var])**(1. / 3)))
    self.xxx = self.dd["x"].reshape((self.n_bins, self.n_bins, self.n_bins))
    self.yyy = self.dd["y"].reshape((self.n_bins, self.n_bins, self.n_bins))
    self.zzz = self.dd["z"].reshape((self.n_bins, self.n_bins, self.n_bins))
    self.center = self.dd.get_field_parameter('center')


  def calculate_projected_cs_eff(self):

    for kk, vv in self.axes.iteritems():

      self.tmp       = self.dd[self.wg_var].reshape((self.n_bins, self.n_bins, self.n_bins))

      self.c_s_eff   = np.sqrt(C_yt.kb * (self.dd['P']+self.dd['P_nt']) / self.dd['density']).reshape((self.n_bins, self.n_bins, self.n_bins))
      self.c_s_eff   = self.c_s_eff.to('km/s')

      # project plane, coordinates.
      self.c_s_eff_proj[kk] = np.sum(self.c_s_eff * self.tmp, axis=int(kk))/np.sum(self.tmp, axis=int(kk))


  def calculate_projected_vel(self):

    """ calculate for all projected planes"""

    self.velocityx = self.dd['velx'].reshape((self.n_bins, self.n_bins, self.n_bins))
    self.velocityy = self.dd['vely'].reshape((self.n_bins, self.n_bins, self.n_bins))
    self.velocityz = self.dd['velz'].reshape((self.n_bins, self.n_bins, self.n_bins))

    for kk, vv in self.axes.iteritems():

      proj = self.ds.proj('velx', int(kk), weight_field=self.wg_var)
      velx_projected = proj['velx']
      self.velx_projected = velx_projected.reshape(
          (-1, int(np.sqrt(velx_projected.shape[0]))))

      proj = self.ds.proj('vely', int(kk), weight_field=self.wg_var)
      vely_projected = proj['vely']
      self.vely_projected = vely_projected.reshape(
          (-1, int(np.sqrt(vely_projected.shape[0]))))

      proj = self.ds.proj('velz', int(kk), weight_field=self.wg_var)
      velz_projected = proj['velz']
      self.velz_projected = velz_projected.reshape(
          (-1, int(np.sqrt(velz_projected.shape[0]))))

      if kk is '2':
          self.vel[kk] = [self.velx_projected, self.vely_projected]
          self.veldisp[kk] = [self.velocityx, self.velocityy]
          self.veldisp_vertical[kk] = np.std(self.velocityz, axis=int(kk))
      elif kk is '1':
          self.vel[kk] = [self.velx_projected, self.velz_projected]
          self.veldisp[kk] = [self.velocityx, self.velocityz]
          self.veldisp_vertical[kk] = np.std(self.velocityy, axis=int(kk))
      elif kk is '0':
          self.vel[kk] = [self.vely_projected, self.velz_projected]
          self.veldisp[kk] = [self.velocityy, self.velocityz]
          self.veldisp_vertical[kk] = np.std(self.velocityx, axis=int(kk))


  def project_SD(self):
    for kk, vv in self.axes.iteritems():

      if self.field_type == 'gas':
        proj = self.ds.proj("density", int(kk), method='integrate')
        self.projected_SurfaceDensity[kk] = proj['density'].reshape(
            (-1, int(np.sqrt(proj['density'].shape[0]))))

      elif self.field_type == 'star':

        # rho of star from mass
        def _rho_star(field, data):
            dx = abs(self.xxx[0][1:] - self.xxx[1][:-1]).min().convert_to_units('pc')
            dy = abs(self.yyy[0][1:] - self.yyy[1][:-1]).min().convert_to_units('pc')
            dz = np.diff(self.zzz[0][0]).min().convert_to_units('pc')

            rho = (data['mass'].convert_to_units('g')) / dx.convert_to_units('cm') / dy.convert_to_units('cm') / dz.convert_to_units('cm')
            return rho

        self.ds.add_field(("rho_star"), function=_rho_star,
              units='g/cm**3')
        _dd = self.ds.all_data()
        if self.megaverbose:
          print _dd['rho_star']
        del _dd

        # mass to surface density
        proj = self.ds.proj('rho_star', int(kk), method='integrate')
        self.projected_SurfaceDensity[kk] = proj['rho_star'].reshape(
            (-1, int(np.sqrt(proj['rho_star'].shape[0]))))


  def project_coords_along_axes(self):
    for kk, vv in self.axes.iteritems():
      if kk is '2':
          self.coords[kk] = [self.xxx[:, :, 0], self.yyy[:, :, 0]]
          self.coords3d[kk] = [self.xxx, self.yyy]
          self.center_plane[kk] = [self.center[0], self.center[1]]
      elif kk is '1':
          self.coords[kk] = [self.xxx[:, 0, :], self.zzz[:, 0, :]]
          self.coords3d[kk] = [self.xxx, self.zzz]
          self.center_plane[kk] = [self.center[0], self.center[2]]
      elif kk is '0':
          self.coords3d[kk] = [self.yyy, self.zzz]
          self.coords[kk] = [self.yyy[0, :, :], self.zzz[0, :, :]]
          self.center_plane[kk] = [self.center[1], self.center[2]]

  def set_smooth_size(self):

    # get spacing in x and y direction in kpc
    nx = np.shape(self.coords_plane[0])[0]
    ny = np.shape(self.coords_plane[1])[1]
    dx = (self.coords_plane[0].max() - self.coords_plane[0].min())/nx
    dy = (self.coords_plane[1].max() - self.coords_plane[1].min())/ny
    dx = (dx.convert_to_units('kpc')).value
    dy = (dy.convert_to_units('kpc')).value
    # set smoothing radius
    self.smooth_size = 0.5*self.smooth_kpc/np.sqrt(dx * dy)

  def project_onto_plane(self):

    """ pick along one plane """

    self.vel_plane              = self.vel[self.plane]
    self.veldisp_plane          = self.veldisp[self.plane]
    self.veldisp_vertical_plane = self.veldisp_vertical[self.plane]
    self.coords_plane           = self.coords[self.plane]
    self.coords3d_plane         = self.coords3d[self.plane]
    self.SD                     = self.projected_SurfaceDensity[self.plane]
    self.set_smooth_size()
    #
    if self.field_type == 'gas':
      self.c_s_eff_plane = self.c_s_eff_proj[self.plane]

    if self.verbose:
      print 'max/min SD', np.max(self.SD), np.min(self.SD)

  def plot_SD(self):

    plt.figure()
    im = plt.imshow(self.SD.value, origin='lower',
                    cmap=cmap)
    plt.title(r'$\Sigma$')
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\Sigma$ [cgs]", fontsize=16)
    if self.show:
      plt.show(block=False)


  def plot_veldisp_vert(self):

    plt.figure()
    im = plt.imshow(self.veldisp_vertical_plane.value, origin='lower',
                    cmap=cmap)
    plt.title('vdisp vertical')
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\sigma_z$ [km s$^{-1}$]", fontsize=16)
    if self.show:
      plt.show(block=False)

  def calc_radial_veldisp(self):
    """ calculate velocity dispersion along radial direction """

    i_hat_3d = self.coords3d_plane[0] - self.center_plane[self.plane][0]
    j_hat_3d = self.coords3d_plane[1] - self.center_plane[self.plane][1]
    R_3d = np.sqrt(i_hat_3d**2 + j_hat_3d**2)
    if self.megaverbose:
      print R_3d.min()
    i_hat_3d /= R_3d
    j_hat_3d /= R_3d

    # weigh by mass
    wg       = self.dd[self.wg_var].reshape((self.n_bins, self.n_bins, self.n_bins))
    v_r      = self.veldisp_plane[0] * i_hat_3d + self.veldisp_plane[1] * j_hat_3d
    mean_v_r = np.sum(wg * v_r, axis=int(self.plane))/np.sum(wg, axis=int(self.plane))
    if self.plane == '0':
      tmp = mean_v_r[np.newaxis,:,:]
    if self.plane == '1':
      tmp = mean_v_r[:,np.newaxis,:]
    if self.plane == '2':
      tmp = mean_v_r[:,:,np.newaxis]

    self.sigma_r  = np.sqrt(np.sum(wg * (v_r - tmp)**2, axis=int(self.plane))/np.sum(wg, axis=int(self.plane)))

    if self.verbose:
      print 'max/min radial vdisp in km/s (from velocity)', np.max(self.sigma_r), np.min(self.sigma_r)

  def plot_radial_veldisp(self):
    plt.figure()
    if self.field_type == 'gas':
      plt.imshow(self.sigma_r, cmap=cmap, origin='lower')
    elif self.field_type == 'star':
      plt.imshow(np.log10(self.sigma_r), cmap=cmap, origin='lower')
    plt.title(r'$\sigma_r$ weighted by '+self.wg_var)
    plt.colorbar()
    if self.show:
      plt.show(block=False)


  def plot_cs_eff_plane(self):
    plt.figure()
    plt.imshow(self.c_s_eff_plane.value, cmap=cmap, origin='lower')
    plt.title(r'$c_{s, {\rm eff}}$ weighted by '+self.wg_var)
    plt.colorbar()
    if self.show:
      plt.show(block=False)
    if self.verbose:
      print r'max/min $c_s$ in km/s (from velocity)', np.max(self.c_s_eff_plane), np.min(self.c_s_eff_plane)


  def calc_tot_sigma_r(self):
    """ apparently Inoue+16 define sigma = sigma_r + pressure contribution (see the text close to eq. 1 in http://adsabs.harvard.edu/abs/2016MNRAS.456.2052I)"""

    if self.field_type == 'gas':
      self.sigma_r = np.sqrt(self.sigma_r**2 + self.c_s_eff_plane**2)
    elif self.field_type == 'star':
      pass      # no changes


  def smooth_sigma_r(self, plot=True):

    # smooth map to regularized the derivates
    if self.smooth_log:
      _sigma_r = 10.**gaussian_filter(np.log10(self.sigma_r), self.smooth_size)
    else:
      _sigma_r = gaussian_filter(self.sigma_r, self.smooth_size)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(121)
      ax.imshow(self.sigma_r.value, cmap=cmap, origin='lower')
      ax = plt.subplot(122)
      ax = plt.subplot(133)
      ax.imshow(_sigma_r, cmap=cmap, origin='lower')
      if self.show:
        plt.show(block=False)

    self.sigma_r = _sigma_r


  def calc_v_phi(self):

    i_hat = self.coords_plane[0] - self.center_plane[self.plane][0]
    j_hat = self.coords_plane[1] - self.center_plane[self.plane][1]
    self.R = np.sqrt(i_hat**2 + j_hat**2)    # + k_hat**2)
    if self.megaverbose:
      print self.R.min()
    self.i_hat = i_hat / self.R
    self.j_hat = j_hat / self.R

    vi = self.vel_plane[0]
    vj = self.vel_plane[1]
    # radial_vel = (vi * i_hat + vj * j_hat)
    # radial_veloDisp = np.std(radial_vel)
    # print 'radial velocity           ', np.max(radial_vel), np.min(radial_vel)
    # print 'radial velocity dispersion', radial_veloDisp

    theta = np.arctan2(self.j_hat, self.i_hat)
    _v_r = np.cos(theta) * vi + np.sin(theta) * vj
    self.v_phi = (-np.sin(theta) * vi + np.cos(theta) * vj)
    if self.verbose:
      print 'maxmin v_phi    ', np.max(self.v_phi), np.min(self.v_phi)
      print 'std v_phi', np.std(self.v_phi)    # km/s


  def smooth_v_phi(self, plot=True):

    if self.smooth_log:
      neg_ind = self.v_phi < 0.0
      _v_phi = 10.**gaussian_filter(np.log10(abs(self.v_phi)), self.smooth_size)
      _v_phi[neg_ind] = - _v_phi[neg_ind]
    else:
      _v_phi = gaussian_filter(self.v_phi, self.smooth_size)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(121)
      ax.imshow(self.v_phi, cmap=cmap, origin='lower')
      ax = plt.subplot(122)
      ax.imshow(_v_phi, cmap=cmap, origin='lower')
      if self.show:
        plt.show(block=False)

    self.v_phi = _v_phi


  def plot_v_phi(self):
    plt.figure()
    if self.field_type == 'gas':
      plt.imshow(self.v_phi, cmap=cmap, origin='lower')
    elif self.field_type == 'star':
      plt.imshow(np.log10(self.v_phi), cmap=cmap, origin='lower')
    plt.colorbar()
    plt.title('vphi')
    if self.show:
      plt.show(block=False)


  def smooth_SD(self, plot=True):
    if self.smooth_log:
      _SD = 10.**gaussian_filter(np.log10(self.SD), self.smooth_size)
    else:
      _SD = gaussian_filter(self.SD, self.smooth_size)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(131)
      ax.imshow(np.log10(self.SD), cmap=cmap, origin='lower')
      ax = plt.subplot(132)
      ax.imshow(np.log10(_SD), cmap=cmap, origin='lower')
      if self.show:
        plt.show(block=False)
    self.SD = _SD


  def calc_omega(self):

    # kappa
    self.R = self.R.value * 1.e3 * self.pc2cm

    try:
      self.v_phi = self.v_phi.value * 1.e5         # cm/s
    except AttributeError:
      self.v_phi = self.v_phi * 1.e5

    cost = self.i_hat
    sint = self.j_hat
    x_slice = (self.coords_plane[0] - self.center_plane[self.plane][0]) * 1.e3 * self.pc2cm
    self.x_slice = x_slice.value
    y_slice = (self.coords_plane[1] - self.center_plane[self.plane][1]) * 1.e3 * self.pc2cm
    self.y_slice = y_slice.value
    self.r_slice = np.sqrt(self.x_slice**2 + self.y_slice**2)    # cm

    self.omega_measured = self.v_phi / self.r_slice
    omega_measured_standard_unit = self.v_phi / 1.e5 / (self.r_slice / self.pc2cm / 1.e3)
    if self.verbose:
      print 'mean omega', (omega_measured_standard_unit).mean(), 'km/s/kpc'      # km/s/kpc
    omega_mw_kms_kpc = 220. / 8
    if self.verbose:
      print 'omega_MW',omega_mw_kms_kpc, 'km/s/kpc'               # ok, unit comparable to MW value

  def plot_omega(self):
    plt.figure()
    plt.imshow(self.omega_measured)
    plt.colorbar()
    if self.show:
      plt.show(block=False)

  def calc_kappa(self, radial_nbins):

    if self.megaverbose:
      print 'calc_kappa()'
    if not radial_nbins:
      radial_nbins = 100

    # Annular bin edges and centers
    bins = np.linspace(0, 1, radial_nbins) * self.r_slice.max()
    bin_centers = bins[:-1] + (bins[1:] - bins[:-1]) / 2.

    # Count how many pixels fall into each radial bin
    hist, _ = np.histogram(self.r_slice, bins)
    hist[hist == 0] = 1
    if self.megaverbose:
      print hist

    # Get flat 1D indices into r_slice for each bin.
    inds = np.digitize(self.r_slice.flat, bins) - 1

    # Calculate mean omega at each radius.
    # Need to append [:-1] to get rid of counts outside bin range.
    omega_of_r = np.bincount(inds, weights=self.omega_measured.flat)[:-1] / hist

    # Calculate standard deviation of the mean omega for each bin.
    omega2_of_r = np.bincount(inds, weights=self.omega_measured.flat[:]**2)[:-1] / hist
    omega_of_r_std = np.sqrt(omega2_of_r - omega_of_r**2)

    # Calculate radial derivative of the rotation frequency using a spline
    # interpolator
    omega = np.zeros(self.r_slice.shape)
    omega_deriv = np.zeros(self.r_slice.shape)

    from scipy import interpolate
    omega_interpolator = interpolate.splrep(bin_centers, omega_of_r, k=5, w=1 / omega_of_r_std)

    omega.flat[:] = interpolate.splev(
        self.r_slice.flat, omega_interpolator)
    rotation_frequency = omega       # should be 1/s

    omega_deriv.flat[:] = interpolate.splev(
        self.r_slice.flat, omega_interpolator, der=1)
    rotation_frequency_derivative = omega_deriv

    domega_dr = np.zeros(self.r_slice.shape)
    domega_dr.flat[:] = interpolate.splev(
        self.r_slice.flat, omega_interpolator, der=1)

    # Finally, calculate epicyclic frequency.
    # plt.figure()
    # plt.imshow(domega_dr)
    # plt.colorbar()
    # plt.show(block=False)

    rotation_frequency = self.omega_measured
    kappa_sq = 2 * rotation_frequency / self.r_slice * (
        2 * self.r_slice * rotation_frequency + self.r_slice**2 * domega_dr)
    kappa_sq[kappa_sq < 0] = np.min(kappa_sq[kappa_sq > 0])
    self.kappa = np.sqrt(kappa_sq)
    if self.megaverbose:
      print self.kappa      # in the MW, kappa ~ omega, which is also true here


  def smooth_kappa(self, plot=True):
    if self.smooth_log:
      _kappa = 10.**gaussian_filter(np.log10(self.kappa), self.smooth_size)
    else:
      _kappa = gaussian_filter(self.kappa, self.smooth_size)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(131)
      ax.imshow(np.log10(self.kappa), cmap=cmap, origin='lower')
      ax = plt.subplot(132)
      ax.imshow(np.log10(_kappa), cmap=cmap, origin='lower')
      if self.show:
        plt.show(block=False)
    self.kappa = _kappa


  def plot_kappa(self):
    plt.figure()
    plt.imshow(np.log10(self.kappa), cmap=cmap, origin='lower')
    plt.colorbar()
    if self.show:
      plt.show(block=False)

  def calc_Q(self):
    # calculate Q_gas
    radial_veloDisp_cgs = self.sigma_r * 1.e5
    whnzero = np.where(self.SD != 0)

    self.Q = np.zeros(self.SD.shape) * np.nan

    if self.field_type == 'gas':
      self.Q[whnzero] = radial_veloDisp_cgs[whnzero] * self.kappa[whnzero] / (self.A_gas * self.G * self.SD[whnzero])

    elif self.field_type == 'star':
      self.Q[whnzero] = radial_veloDisp_cgs[whnzero] * self.kappa[whnzero] / \
                        (self.A_star * self.G * self.SD[whnzero])

    if self.megaverbose:
      print "Q of field {0:s}: {1:}".format(self.field_type, self.Q)
    # print(np.isnan(self.Q) == True).any()

    return self.Q


  def set_plot_range_set_by_camera(self):
    self._xmin = (self.coords_plane[0] - self.center_plane[self.plane][0]).min()
    self._xmax = (self.coords_plane[0] - self.center_plane[self.plane][0]).max()
    self._ymin = (self.coords_plane[1] - self.center_plane[self.plane][1]).min()
    self._ymax = (self.coords_plane[1] - self.center_plane[self.plane][1]).max()


  def plot_Q(self):
    plt.figure()
    plt.imshow(np.log10(self.Q), origin='lower',
                extent=(self._xmin, self._xmax, self._ymin, self._ymax)
                #, cmap=cmap
                )
    plt.title('Log Q without blanking for field: {:s}'.format(self.field_type))
    cbar = plt.colorbar()
    cbar.set_label(r"$\log{Q}$", fontsize=16)
    if self.show:
      plt.show(block=False)


  def plot_all_quant(self):

    # show all quantities in one figure
    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(left=0.10, right=0.90, hspace=0.1, wspace=0.25)
    ax = plt.subplot(221)
    im = ax.imshow(np.log10(self.SD / self.cm2pc**2 * self.g2Msun),
                    origin='lower', extent=(self._xmin, self._xmax, self._ymin, self._ymax),
                    cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]", fontsize=16)

    ax = plt.subplot(222)
    im = ax.imshow(self.sigma_r, origin='lower', extent=(self._xmin, self._xmax, self._ymin, self._ymax), cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\sigma$ [km\,s$^{-1}$]", fontsize=16)
    if self.show:
      plt.show(block=False)

    ax = plt.subplot(223)
    im = ax.imshow(np.log10(self.kappa * 3.086e+16), origin='lower', extent=(self._xmin, self._xmax, self._ymin, self._ymax), cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]", fontsize=16)

    ax = plt.subplot(224)
    import matplotlib as mpl
    im = ax.imshow(np.log10(self.Q), origin='lower', \
                   extent=(self._xmin, self._xmax, self._ymin, self._ymax),
                   cmap=cmap_div,
                   vmin=-1,
                   vmax= 1
                   )
    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
    # cbar.set_clim(-1, 1)                      # set color max min range
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q}$", fontsize=16)

    # plt.tight_layout()
    if self.show:
      plt.show(block=False)


  def plot_all_quant_zoom(self, central_kpc_one_side=None, annotate_clump=True,
                          clump_list_filename=None):
    """
      show all quantities in one figure, but only show central region of the plot
      (i.e., on the main galaxy)
      from -1.5 kpc to 1.5 kpc
    """

    if not central_kpc_one_side:
      central_kpc_one_side = 1.5

    if annotate_clump:
      assert clump_list_filename is not None

    xspacing = (self._xmax - self._xmin)/len(self.Q)
    xruler = np.arange(self._xmin, self._xmax, xspacing)
    rightBound = np.argmin(abs(xruler - central_kpc_one_side))
    leftBound = np.argmin(abs(xruler + central_kpc_one_side))

    yspacing = (self._ymax - self._ymin)/len(self.Q)
    yruler = np.arange(self._ymin, self._ymax, yspacing)
    topBound = np.argmin(abs(yruler - central_kpc_one_side))
    bottomBound = np.argmin(abs(yruler + central_kpc_one_side))

    x1, x2 = xruler[leftBound],  xruler[rightBound]
    y1, y2 = yruler[bottomBound],yruler[topBound]

    fig = plt.figure(figsize=(9, 8))
    fig.subplots_adjust(left=0.10, right=0.90, hspace=0.3, wspace=0.25)
    ax = plt.subplot(221)
    im = ax.imshow(np.log10(self.SD / self.cm2pc**2 * self.g2Msun)[bottomBound: topBound, leftBound:rightBound],
              origin='lower',
              extent=(xruler[leftBound],
                      xruler[rightBound],
                      yruler[bottomBound],
                      yruler[topBound]),
              cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]", fontsize=16)

    if annotate_clump:
      _, clx, cly, clz = np.loadtxt(clump_list_filename, unpack=True)
      # pos2 = pos * self.factor_R

      if self.plane == '0':
        plt.plot(clz, cly, '*', markersize=15, markerfacecolor='none',
                   markeredgecolor=self.c_clump, fillstyle='none',
                   markeredgewidth=1.15,
                  linewidth=1.8)
      elif self.plane == '1':
        plt.plot(clz, clx, '*', markersize=15, markerfacecolor='none',
                   markeredgecolor=self.c_clump, fillstyle='none',
                   markeredgewidth=1.15,
                  linewidth=1.8)
      elif self.plane == '2':
        plt.plot(cly, clx, '*', markersize=15, markerfacecolor='none',
                   markeredgecolor=self.c_clump, fillstyle='none',
                   markeredgewidth=1.15,
                  linewidth=1.8)
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)

    ax = plt.subplot(222)
    im = ax.imshow(self.sigma_r[bottomBound: topBound, leftBound:rightBound],
                   origin='lower',
                   extent=(xruler[leftBound],
                           xruler[rightBound],
                           yruler[bottomBound],
                           yruler[topBound]),
                   cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\sigma$ [km\,s$^{-1}$]", fontsize=16)
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)

    ax = plt.subplot(223)
    im = ax.imshow(np.log10(self.kappa * 3.086e+16)[bottomBound: topBound, leftBound:rightBound],
                   origin='lower',
                   extent=(xruler[leftBound],
                           xruler[rightBound],
                           yruler[bottomBound],
                           yruler[topBound]),
                   cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]", fontsize=16)
    plt.xlabel('kpc', fontsize=16)
    plt.ylabel('kpc', fontsize=16)
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)

    ax = plt.subplot(224)
    map_Q = np.log10(gaussian_filter(self.Q, sigma=self.smooth_size))
    im = ax.imshow(map_Q[bottomBound: topBound, leftBound:rightBound],
                   origin='lower',
                   extent=(xruler[leftBound],
                           xruler[rightBound],
                           yruler[bottomBound],
                           yruler[topBound]),
                   cmap=cmap_div,
                   vmin=-1, vmax=1     # clip at -1 < log10(Q) < 1
                   )
    if annotate_clump:
      _, clx, cly, clz = np.loadtxt(clump_list_filename, unpack=True)
      # pos2 = pos * self.factor_R

      if self.plane == '0':
        plt.plot(clz, cly, '*', markersize=15, markerfacecolor='none',
                   markeredgecolor=self.c_clump, fillstyle='none',
                   markeredgewidth=1.65,
                  linewidth=1.8)
      elif self.plane == '1':
        plt.plot(clz, clx, '*', markersize=15, markerfacecolor='none',
                   markeredgecolor=self.c_clump, fillstyle='none',
                   markeredgewidth=1.65,
                  linewidth=1.8)
      elif self.plane == '2':
        plt.plot(cly, clx, '*', markersize=15, markerfacecolor='none',
                   markeredgecolor=self.c_clump, fillstyle='none',
                   markeredgewidth=1.65,
                  linewidth=1.8)

    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q}$", fontsize=16)

    plt.tight_layout()
    if self.show:
      plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_' + self.field_type + '_toomre_proj_' + self.plane +\
            '_zoom_'+str(central_kpc_one_side)+'_kpc'+\
            '.png'
    if self.verbose:
      print 'save to'
      print '  ',self.fold_out+out_f
    plt.savefig(self.fold_out+out_f)
    plt.clf()


  def plot_single_SD(self, central_kpc_one_side=None):
    """
    plot and save single panel SD

    """

    if not central_kpc_one_side:
      central_kpc_one_side = 1.5

    xspacing = (self._xmax - self._xmin)/len(self.SD)
    xruler = np.arange(self._xmin, self._xmax, xspacing)
    rightBound = np.argmin(abs(xruler - central_kpc_one_side))
    leftBound = np.argmin(abs(xruler + central_kpc_one_side))

    yspacing = (self._ymax - self._ymin)/len(self.SD)
    yruler = np.arange(self._ymin, self._ymax, yspacing)
    topBound = np.argmin(abs(yruler - central_kpc_one_side))
    bottomBound = np.argmin(abs(yruler + central_kpc_one_side))

    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(left=0.10, right=0.90, bottom=0.15, top=0.9)
    ax = plt.gca()
    im = ax.imshow(np.log10(self.SD / self.cm2pc**2 * \
                            self.g2Msun)[bottomBound: topBound,
                                         leftBound:rightBound],
                   origin='lower',
                   extent=(xruler[leftBound],
                           xruler[rightBound],
                           yruler[bottomBound],
                           yruler[topBound]),
                   cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(r"$\log{\Sigma}$ [M$_{\odot}$~pc$^{-2}$]", fontsize=16)

    plt.tight_layout()
    if self.show:
      plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_' + self.field_type + '_SD_proj_' + \
            self.plane + '_zoom_' + str(central_kpc_one_side) + '_kpc' + '.png'

    if self.verbose:
      print 'save to'
      print '  ',self.fold_out+out_f
    plt.savefig(self.fold_out+out_f, bbox_inches='tight')
    plt.clf()


  def run(self, radial_nbins=None, central_kpc_one_side=None,
          annotate_clump=False, clump_list_filename=None):
    self.load_cam_stuff()
    self.load_data()
    self.setup_config()

    if self.field_type == 'gas':
      self.calculate_projected_cs_eff()

    self.calculate_projected_vel()
    self.project_SD()
    self.project_coords_along_axes()
    self.project_onto_plane()
    self.calc_radial_veldisp()
    self.calc_tot_sigma_r()
    self.smooth_sigma_r(plot=self.debug)
    self.calc_v_phi()
    self.smooth_v_phi(plot=self.debug)
    self.smooth_SD(plot=self.debug)
    self.calc_omega()
    self.calc_kappa(radial_nbins)
    self.smooth_kappa(plot=self.debug)
    self.calc_Q()

    self.set_plot_range_set_by_camera()
    self.plot_all_quant_zoom(central_kpc_one_side, annotate_clump,
                             clump_list_filename)

    self.plot_single_SD(central_kpc_one_side)

    if self.debug:
      self.plot_SD()
      self.plot_veldisp_vert()
      self.plot_radial_veldisp()
      if self.field_type == 'gas':
        self.plot_cs_eff_plane()
      self.plot_v_phi()
      self.plot_omega()
      self.plot_kappa()
      self.plot_Q()

    return self.Q



class ToomreAnalyze_2comp(object):
  """

  Toomre for gas + star in thick disk

  """
  def __init__(self, Q_gas, Q_star,fold_out='' ,smooth_size = 0):

    """

    Calculates the combined Toomre Q parameter

    Take into account the finite thickness of the disk
    and that the gas and stars independently contribute to the gravitational
    potential.

    See Romeo & Wiegert (2011) [2011MNRAS.416.1191R] for details
    http://adsabs.harvard.edu/abs/2011MNRAS.416.1191R

    Parameters
    ----------
    Q_gas and Q_star are objects


    """

    # regulate output
    self.fold_out       = fold_out
    if (self.fold_out !='' and not os.path.isdir(self.fold_out)):
      os.mkdir(self.fold_out)

    self.debug       = Q_gas.debug or Q_star.debug
    self.verbose     = Q_gas.verbose or Q_star.verbose
    self.megaverbose = Q_gas.megaverbose or Q_star.megaverbose
    self.show        = Q_gas.show or Q_star.show

    self.Q_gas       = Q_gas
    self.Q_star      = Q_star

    assert self.Q_gas.plane == self.Q_star.plane
    self.plane = self.Q_gas.plane

    assert self.Q_gas.isnap == self.Q_star.isnap
    self.isnap = self.Q_gas.isnap

    self.smooth_size = smooth_size
    self.c_clump = self.Q_star.c_clump

  def compute_T(self, veldisp_vert, veldisp_r):
    """ see eq. 4 of Inoue+16"""


    ratio       = np.ones_like(veldisp_vert)
    mask        = veldisp_r > 0
    ratio[mask] = veldisp_vert[mask]/veldisp_r[mask]

    if self.verbose:
      print "vel disp ratio (sigma_z/sigma_r): "
      print '  max/ min',np.max(ratio),np.min(ratio)
      print '  mean/std',np.mean(ratio),np.std(ratio)

    res1 = 1.  + 0.6 * ratio**2
    res2 = 0.8 + 0.7 * ratio

    out            = res1
    out[ratio>0.5] = res2[ratio>0.5]

    return out

  def compute_T_s(self):
    self.T_s = self.compute_T(self.Q_star.veldisp_vertical_plane.value, self.Q_star.sigma_r)


  def compute_T_g(self):
    self.T_g = self.compute_T(self.Q_gas.veldisp_vertical_plane.value, self.Q_gas.sigma_r)


  def calc_Q_eff(self):

    # The effect of thickness is to increase the stability parameter of each
    # component by a factor T, which depends on the ratio of vertical to
    # radial velocity dispersion.

    """ see eq. 3 of Inoue+16"""
    # calculate weight
    up      = 2. * self.Q_star.sigma_r * self.Q_gas.sigma_r
    low     = (self.Q_star.sigma_r**2 + self.Q_gas.sigma_r**2)
    w       = np.zeros_like(up)
    mask    = low>0
    w[mask] = up[mask]/low[mask]

    # get thick Q parameters
    Q_thick_star = self.Q_star.Q * self.T_s
    Q_thick_gas  = self.Q_gas.Q  * self.T_g

    # regularize
    Q_thick_star[np.isnan(Q_thick_star)] = 0
    Q_thick_gas [np.isnan(Q_thick_gas )] = 0

    if self.verbose:
      for arr,nome in zip([w,Q_thick_star, Q_thick_gas],['weight','Q_thick_star','Q_thick_gas']):
        print nome
        print '  max/ min',np.max(arr),np.min(arr)
        print '  mean/std',np.mean(arr),np.std(arr)

    # init
    mask          =  np.logical_and(Q_thick_star> 0, Q_thick_gas > 0)
    Q_twoComp_inv = - np.ones_like(Q_thick_gas)
    self.Q_twoComp= np.zeros_like(Q_thick_gas)
    res1,res2     = - np.ones_like(Q_thick_gas), - np.ones_like(Q_thick_gas)

    # compute Q composite
    """ see eq. 3 of Inoue+16"""
    res1[mask]    = w[mask] /Q_thick_star[mask] + 1.     / Q_thick_gas[mask]
    res2[mask]    = 1.      /Q_thick_star[mask] + w[mask]/Q_thick_gas[mask]
    #
    mask2                = np.logical_and(Q_thick_star > Q_thick_gas,mask)
    Q_twoComp_inv[mask2] = res1[mask2]
    mask2                = np.logical_and(Q_thick_star < Q_thick_gas,mask)
    Q_twoComp_inv[mask2] = res1[mask2]
    #
    mask                 =     Q_twoComp_inv > 0
    self.Q_twoComp[mask] = 1. / Q_twoComp_inv[mask]

    if self.verbose:

      for arr,nome in zip([self.Q_twoComp],['Q_2_comp']):
        print nome
        print '  max/ min',np.max(arr),np.min(arr)
        print '  mean/std',np.mean(arr),np.std(arr)

      if (np.isnan(self.Q_twoComp) == True).any():
        print 'NANs are present'

    return self.Q_twoComp


  def smooth_Q_eff(self):
    self.smooth_size = self.Q_gas.smooth_size
    self.Q_twoComp = gaussian_filter(self.Q_twoComp, sigma=self.smooth_size)


  def plot_Q_eff(self):
    plt.figure(figsize=(5, 5))
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.9)
    ax = plt.gca()
    im = plt.imshow(np.log10(self.Q_twoComp),
                    extent=(self.Q_gas._xmin, self.Q_gas._xmax, self.Q_gas._ymin, self.Q_gas._ymax),
                    cmap=cmap_div,
                    origin='lower',
                    vmin=-1, vmax=1)     # clip at -1 < log10(Q) < 1
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1],
                         cax=cax
                        )
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q_{\rm eff}}$", fontsize=16)
    plt.title(r'$Q_{\rm eff}$')
    plt.xlabel('kpc', fontsize=16)
    plt.ylabel('kpc', fontsize=16)

    plt.tight_layout()
    if self.show:
      plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_toomreEff_proj_' + self.plane + '.png'
    if self.verbose:
      print 'save to'
      print '  ',self.fold_out+out_f
    plt.savefig(self.fold_out+out_f, bbox_inches='tight')
    plt.clf()


  def get_bound_zoom(self,central_kpc_one_side = None):

    if central_kpc_one_side is None:
      central_kpc_one_side = 1.5

    xspacing = (self.Q_gas._xmax - self.Q_gas._xmin)/len(self.Q_twoComp)
    xruler = np.arange(self.Q_gas._xmin, self.Q_gas._xmax, xspacing)
    rightBound = np.argmin(abs(xruler - central_kpc_one_side))
    leftBound = np.argmin(abs(xruler + central_kpc_one_side))

    yspacing = (self.Q_gas._ymax - self.Q_gas._ymin)/len(self.Q_twoComp)
    yruler = np.arange(self.Q_gas._ymin, self.Q_gas._ymax, yspacing)
    topBound = np.argmin(abs(yruler - central_kpc_one_side))
    bottomBound = np.argmin(abs(yruler + central_kpc_one_side))

    return bottomBound, topBound, leftBound,rightBound,xruler,yruler 

  def plot_Q_eff_zoom(self, central_kpc_one_side=None, annotate_clump=False, clump_list_filename=None):

    if annotate_clump:
      assert clump_list_filename is not None

    bottomBound, topBound, leftBound,rightBound,xruler,yruler = self.get_bound_zoom(central_kpc_one_side= central_kpc_one_side)

    plt.figure()
    # fig.subplots_adjust(left=0.10, right=0.90, hspace=0.3, wspace=0.25)
    im = plt.imshow(np.log10(self.Q_twoComp)[bottomBound: topBound,
                                             leftBound:rightBound],
                    origin='lower',
                    extent=(xruler[leftBound],
                            xruler[rightBound],
                            yruler[bottomBound],
                            yruler[topBound]),
                    cmap=cmap_div,
                    vmin=-1, vmax=1)     # clip at -1 < log10(Q) < 1
    if annotate_clump:
      _, clx, cly, clz = np.loadtxt(clump_list_filename, unpack=True)
      # pos2 = pos * self.factor_R

      if self.plane == '0':
        plt.plot(clz, cly, '*', markersize=15, markerfacecolor='none',
                 markeredgecolor=self.c_clump, fillstyle='none',
                 markeredgewidth=1.65,
                 linewidth=1.8)
      elif self.plane == '1':
        plt.plot(clz, clx, '*', markersize=15, markerfacecolor='none',
                 markeredgecolor=self.c_clump, fillstyle='none',
                 markeredgewidth=1.65,
                 linewidth=1.8)
      elif self.plane == '2':
        plt.plot(cly, clx, '*', markersize=15, markerfacecolor='none',
                 markeredgecolor=self.c_clump, fillstyle='none',
                 markeredgewidth=1.65,
                 linewidth=1.8)

    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q_{\rm eff}}$", fontsize=16)
    plt.title(r'$Q_{\rm eff}$')
    plt.xlabel('kpc', fontsize=16)
    plt.ylabel('kpc', fontsize=16)

    plt.tight_layout()
    if self.show:
      plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_toomreEff_proj_' + self.plane + \
            '_zoom_'+ str(central_kpc_one_side) + '_kpc.png'
    if self.verbose:
      print 'save to'
      print '  ',self.fold_out+out_f
    plt.savefig(self.fold_out+out_f, bbox_inches='tight')
    plt.clf()

  def smooth_maps(self):

    self.Q_star.sigma_r                     = gaussian_filter(self.Q_star.sigma_r                     ,self.smooth_size)
    self.Q_gas.sigma_r                      = gaussian_filter(self.Q_gas.sigma_r                      ,self.smooth_size)

    from yt import YTArray

    tmp = gaussian_filter(self.Q_gas.veldisp_vertical_plane.convert_to_units('km/s').value ,self.smooth_size)
    self.Q_gas.veldisp_vertical_plane = YTArray(tmp, 'km/s')

    tmp = gaussian_filter(self.Q_star.veldisp_vertical_plane.convert_to_units('km/s').value,self.smooth_size)
    self.Q_star.veldisp_vertical_plane = YTArray(tmp, 'km/s')

  def wrapper_plot(self,ax= None
                  ,mappa= None, label= '', tipo= ''
                  ,bottomBound=0, topBound=-1, leftBound=0,rightBound=-1
                  ,xruler= None,yruler= None 
                  ):


    cmap_plt = cmap
    vmin_plt = None
    vmax_plt = None

    if tipo == 'Q':
      cmap_plt=cmap_div
      vmin_plt=-1
      vmax_plt =1
    elif tipo == 'Sigma':
      cmap_plt = cmap2

    im = ax.imshow(mappa[bottomBound: topBound, leftBound:rightBound],
              origin='lower',
              extent=(xruler[leftBound],xruler[rightBound],yruler[bottomBound],yruler[topBound]),
              cmap=cmap_plt, vmin= vmin_plt,vmax=vmax_plt)
    if tipo == 'Q':
      cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
      cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])

    else:
      cbar = plt.colorbar(im)
    cbar.set_label(label, fontsize=16)

    return ax,cbar

  def wrapper_annotate_clumps(self, ax= None, clump_list_filename= None):

    __, clx, cly, clz = np.loadtxt(clump_list_filename, unpack=True)
    # pos2 = pos * self.factor_R
    if self.plane == '0':
      ax.plot(clz, cly, '*', markersize=15, markerfacecolor='none',
                 markeredgecolor=self.Q_gas.c_clump, fillstyle='none',
                 markeredgewidth=1.15,
                linewidth=1.8)
    elif self.plane == '1':
      ax.plot(clz, clx, '*', markersize=15, markerfacecolor='none',
                 markeredgecolor=self.Q_gas.c_clump, fillstyle='none',
                 markeredgewidth=1.15,
                linewidth=1.8)
    elif self.plane == '2':
      ax.plot(cly, clx, '*', markersize=15, markerfacecolor='none',
                 markeredgecolor=self.Q_gas.c_clump, fillstyle='none',
                 markeredgewidth=1.15,
                linewidth=1.8)
    return ax

  def plots_combined(self,central_kpc_one_side = 1.5, annotate_clump = False, clump_list_filename = None,
                     type_plots = '2by2'):

    if not central_kpc_one_side:
      central_kpc_one_side = 1.5

    bottomBound, topBound, leftBound,rightBound,xruler,yruler = self.get_bound_zoom(central_kpc_one_side= central_kpc_one_side)
    x1, x2      = xruler[leftBound],  xruler[rightBound]
    y1, y2      = yruler[bottomBound],yruler[topBound]


    if type_plots in ['2by2','3by2']:

      fig = plt.figure(figsize=(9, 8))
      fig.subplots_adjust(left=0.10, right=0.90, hspace=0.3, wspace=0.25)

      list_maps    = [
                      np.log10(self.Q_gas.SD / self.Q_gas.cm2pc**2 * self.Q_gas.g2Msun)
                     ,np.log10(self.Q_star.SD / self.Q_star.cm2pc**2 * self.Q_star.g2Msun)
                     ,np.log10(self.Q_gas.sigma_r)   
                     ,np.log10(self.Q_star.sigma_r)   
                     ]
      list_plt_ids = [221,222,223,224]
      list_types   = ['Sigma','Sigma','sigma','sigma']
      list_labels  = [
                      r"$\log{\Sigma_{\rm gas}}$ [M$_{\odot}$~pc$^{-2}$]"
                     ,r"$\log{\Sigma_\star}$ [M$_{\odot}$~pc$^{-2}$]"
                     ,r"$\log{\sigma_{\rm gas}}$ [${\rm km}\,{\rm s}^{-1}$]"
                     ,r"$\log{\sigma_\star}$ [${\rm km}\,{\rm s}^{-1}$]"
                     ]
      if type_plots == '3by2':
        list_plt_ids = [321,322,323,324,325,326]
        list_types   = list_types + ['Q','Q']
        list_maps    = list_maps  + [np.log10(self.Q_gas.Q),np.log10(self.Q_star.Q)] 
        list_labels  = list_labels+ [r'$\log Q_{\rm gas}$','$\log Q_{\star}$']
    elif type_plots == '3by1':

      fig = plt.figure(figsize=(9, 3))
      fig.subplots_adjust(left=0.10, right=0.90, hspace=0.3, wspace=0.25)

      list_plt_ids = [131,132,133]
      list_types   = ['Q','Q','Q']
      list_maps    = [np.log10(self.Q_gas.Q),np.log10(self.Q_star.Q),np.log10(self.Q_twoComp)] 
      list_labels  = [r'$\log Q_{\rm gas}$','$\log Q_{\star}$',r'$\log Q_{\rm eff}$']

    for mappa, label, tipo, plt_id in zip(list_maps,list_labels,list_types,list_plt_ids):

      ax     = plt.subplot(plt_id)
      ax, cb = self.wrapper_plot(ax=ax,mappa = mappa,label = label,tipo  = tipo
                           ,bottomBound=bottomBound, topBound=topBound, leftBound= leftBound,rightBound = rightBound
                           ,xruler = xruler , yruler = yruler 
                           )
      if annotate_clump:
        assert clump_list_filename is not None
        ax= self.wrapper_annotate_clumps(ax=ax, clump_list_filename=clump_list_filename)
      ax.set_xlim(x1,x2)
      ax.set_ylim(y1,y2)

    '''
    np.log10(self.Q_gas.kappa * 3.086e+16)
    r"$\log{\kappa}$ [km\,s$^{-1}$\,kpc$^{-1}$]"
    self.Q_gas.Q
    '''

    plt.tight_layout()
    if self.show:
      plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_toomre_combined_'+type_plots+'_' + self.plane +\
            '.png'
    if self.verbose:
      print 'save to'
      print '  ',self.fold_out+out_f
    plt.savefig(self.fold_out+out_f)
    plt.clf()


  def run(self, central_kpc_one_side, annotate_clump, clump_list_filename):

    if self.smooth_size > 0:
      self.smooth_maps()

    self.compute_T_g()
    self.compute_T_s()
    self.calc_Q_eff()
    self.smooth_Q_eff()

    self.plot_Q_eff_zoom(central_kpc_one_side, annotate_clump, clump_list_filename)


if __name__ == '__main__':

  plt.close('all')

  base_out = 'out_toomre/'

  if not os.path.isdir(base_out):
    os.mkdir(base_out)

  plane     = '0'
  isnap     = 16
  annotate  = True

  smooth_kpc  = 0.05   # 50 pc
  clump_cut   = 6.81

  min_mass  = 1.e+1 # used to clip 0 in the stellar mass field
  size_kpc  = 2.0

  #testfile  = 'ss'+str(isnap)+'_h2density_clumppos_ncut_'+str(clump_cut)+'_Ncellmin_10.txt'
  testfile  = 'ss'+str(isnap)+'_h2density_wgclumppos_ncut_'+str(clump_cut)+'_Ncellmin_10.txt'
  fold_out  = base_out+'snap_'+str(isnap)+'/'


  Q_gas_obj = ToomreAnalyze(isnap=isnap, wg_var='density',field_type='gas', plane=plane
                            ,smooth_size_kpc =smooth_kpc
                            ,debug=False
                            ,fold_out = fold_out
                            )

  Q_gas_val = Q_gas_obj.run(radial_nbins=100, central_kpc_one_side=size_kpc
                           ,annotate_clump=annotate,clump_list_filename=testfile
                           )
  Q_gas_obj.plot_all_quant_zoom(size_kpc, annotate_clump=annotate,clump_list_filename=testfile)
  Q_gas_obj.plot_single_SD(central_kpc_one_side=size_kpc)

  Q_star_obj = ToomreAnalyze(isnap=isnap, wg_var='mass',field_type='star', plane=plane
                       ,smooth_size_kpc=smooth_kpc,min_wg = min_mass
                       ,debug=False
                       ,fold_out = fold_out
                        )
  Q_star_val = Q_star_obj.run(radial_nbins=100, central_kpc_one_side=size_kpc
                             ,annotate_clump=annotate,clump_list_filename=testfile
                             )
  Q_star_obj.plot_single_SD(central_kpc_one_side=size_kpc)

  Q_tot_obj = ToomreAnalyze_2comp(Q_gas_obj, Q_star_obj,fold_out = fold_out, smooth_size =Q_gas_obj.smooth_size )
  Q_tot_val = Q_tot_obj.run(central_kpc_one_side=size_kpc
                           ,annotate_clump=annotate,
                           clump_list_filename=testfile
                           )

  for type_plots in ['2by2', '3by2', '3by1']:
    Q_tot_obj.plots_combined(
                            central_kpc_one_side = size_kpc
                            ,annotate_clump=annotate
                            ,clump_list_filename=testfile
                            ,type_plots= type_plots
                           )

  #Q_tot_obj.plot_Q_eff_zoom(central_kpc_one_side=size_kcp,
  #                          annotate_clump=annotate,
  #                          clump_list_filename=testfile)



