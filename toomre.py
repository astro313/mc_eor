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
cmap = cm.get_cmap('viridis')          # 'magma'
cmap_div = cm.get_cmap('RdBu')         # divergent cmap

from pymses.utils import constants as C_py
from yt import units as C_yt
import astropy.constants as C_ap
import astropy.units as U


class ToomreAnalyze(object):
  """
  Single Toomre Q object

  """
  def __init__(self, isnap, wg_var, field_type, plane, bin_size_in_log10=0.1, read_proper_unit=True, verbose=True, debug=False, convertPart=True, megaverbose = False):


    self.isnap = isnap
    self.read_proper_unit = read_proper_unit
    self.wg_var = wg_var         # density for gas, mass for stars
    self.debug = debug
    self.verbose = verbose
    self.megaverbose = megaverbose
    self.convertPart = convertPart
    self.plane = plane
    self.field_type = field_type
    #
    self.bin_size_in_log10 = bin_size_in_log10     # bin to smooth in log10 space
    self.smooth_log = False # i do think we should smooth in linear space

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
    # axis : corresponding to the axis to slice along

    self.axes = {'0': 'x', '1': 'y', '2': 'z'}
    self.vel = {}
    self.veldisp = {}
    self.veldisp_vertical = {}
    self.c_s_eff_proj = {}
    self.projected_SurfaceDensity = {}
    self.coords = {}
    self.coords3d = {}
    self.center_plane = {}


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
      self.data = import_fetch_gal(isnap=self.isnap)
      if self.megaverbose: 
        print self.data.keys()
      self.ds, self.dd = prepare_unigrid(data=self.data,
                               add_unit=True,
                               regionsize_kpc=self.region_size_kpc,
                               debug=self.debug)

    elif self.field_type == 'star':
      # read from resampled.h5
      self.starData = import_fetch_stars(isnap=self.isnap,
                                         verbose=self.verbose,
                                         convert=self.convertPart)
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

    for kk, vv in self.axes.iteritems():

      self.velocityx = self.dd['velx'].reshape((self.n_bins, self.n_bins, self.n_bins))
      self.velocityy = self.dd['vely'].reshape((self.n_bins, self.n_bins, self.n_bins))
      self.velocityz = self.dd['velz'].reshape((self.n_bins, self.n_bins, self.n_bins))

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


  def project_onto_plane(self):

    """ pick along one plane """

    self.vel_plane = self.vel[self.plane]
    self.veldisp_plane = self.veldisp[self.plane]
    self.veldisp_vertical_plane = self.veldisp_vertical[self.plane]
    self.coords_plane = self.coords[self.plane]
    self.coords3d_plane = self.coords3d[self.plane]
    self.SD = self.projected_SurfaceDensity[self.plane]

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
    plt.show(block=False)


  def plot_veldisp_vert(self):

    plt.figure()
    im = plt.imshow(self.veldisp_vertical_plane.value, origin='lower',
                    cmap=cmap)
    plt.title('vdisp vertical')
    cbar = plt.colorbar(im)
    cbar.set_label(r"$\sigma_z$ [km s$^{-1}$]", fontsize=16)
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
    plt.show(block=False)


  def plot_cs_eff_plane(self):
    plt.figure()
    plt.imshow(self.c_s_eff_plane.value, cmap=cmap, origin='lower')
    plt.title(r'$c_{s, {\rm eff}}$ weighted by '+self.wg_var)
    plt.colorbar()
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
      _sigma_r = 10.**gaussian_filter(np.log10(self.sigma_r), self.bin_size_in_log10)
    else:
      _sigma_r = gaussian_filter(self.sigma_r, self.bin_size_in_log10)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(121)
      ax.imshow(self.sigma_r.value, cmap=cmap, origin='lower')
      ax = plt.subplot(122)
      ax = plt.subplot(133)
      ax.imshow(_sigma_r, cmap=cmap, origin='lower')
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
      _v_phi = 10.**gaussian_filter(np.log10(abs(self.v_phi)), self.bin_size_in_log10)
      _v_phi[neg_ind] = - _v_phi[neg_ind]
    else:
      _v_phi = gaussian_filter(self.v_phi, self.bin_size_in_log10)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(121)
      ax.imshow(self.v_phi, cmap=cmap, origin='lower')
      ax = plt.subplot(122)
      ax.imshow(_v_phi, cmap=cmap, origin='lower')
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
    plt.show(block=False)


  def smooth_SD(self, plot=True):
    if self.smooth_log:
      _SD = 10.**gaussian_filter(np.log10(self.SD), self.bin_size_in_log10)
    else:
      _SD = gaussian_filter(self.SD, self.bin_size_in_log10)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(131)
      ax.imshow(np.log10(self.SD), cmap=cmap, origin='lower')
      ax = plt.subplot(132)
      ax.imshow(np.log10(_SD), cmap=cmap, origin='lower')
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
      _kappa = 10.**gaussian_filter(np.log10(self.kappa), self.bin_size_in_log10)
    else:
      _kappa = gaussian_filter(self.kappa, self.bin_size_in_log10)

    if plot:
      fig = plt.figure()
      ax = plt.subplot(131)
      ax.imshow(np.log10(self.kappa), cmap=cmap, origin='lower')
      ax = plt.subplot(132)
      ax.imshow(np.log10(_kappa), cmap=cmap, origin='lower')
      plt.show(block=False)
    self.kappa = _kappa


  def plot_kappa(self):
    plt.figure()
    plt.imshow(np.log10(self.kappa), cmap=cmap, origin='lower')
    plt.colorbar()
    plt.show(block=False)


  def calc_Q(self):
    # calculate Q_gas
    radial_veloDisp_cgs = self.sigma_r * 1.e5
    whnzero = np.where(self.SD != 0)

    self.Q = np.zeros(self.SD.shape) * np.nan

    if self.field_type == 'gas':
      self.Q[whnzero] = radial_veloDisp_cgs[whnzero] * self.kappa[whnzero] / (self.A_gas * self.G * self.SD[whnzero])

    elif self.field_type == 'star':
      self.Q[whnzero] = radial_veloDisp_cgs * self.kappa[whnzero] / \
                        (self.A_star * self.G * self.SD[whnzero])

    if self.megaverbose: 
      print "Q of field {0:s}: {1:}".format(self.field_type, self.Q)
    # print(np.isnan(self.Q) == True).any()

    return self.Q


  def plot_range_set_by_camera(self):
    self._xmin = (self.coords_plane[0] - self.center_plane[self.plane][0]).min()
    self._xmax = (self.coords_plane[0] - self.center_plane[self.plane][0]).max()
    self._ymin = (self.coords_plane[1] - self.center_plane[self.plane][1]).min()
    self._ymax = (self.coords_plane[1] - self.center_plane[self.plane][1]).max()


  def plot_Q(self):
    # define plot boundary kpc

    plt.figure()
    plt.imshow(np.log10(self.Q), origin='lower',
                extent=(self._xmin, self._xmax, self._ymin, self._ymax)
                #, cmap=cmap
                )
    plt.title('Log Q without blanking for field: {:s}'.format(self.field_type))
    cbar = plt.colorbar()
    cbar.set_label(r"$\log{Q}$", fontsize=16)
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
    plt.show(block=False)

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

    ax = plt.subplot(224)

    map_Q = np.log10(gaussian_filter(self.Q, sigma=0.55))

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
        plt.plot(cly, clz, 'x', markersize=15, color='darkturquoise')
      elif self.plane == '1':
        plt.plot(clx, clz, 'x', markersize=15, color='antiquewhite')
      elif self.plane == '2':
        plt.plot(clx, cly, 'x', markersize=15, color='antiquewhite')

    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q}$", fontsize=16)

    # plt.tight_layout()
    plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_' + self.field_type + '_toomre_proj_' + self.plane +\
            '_zoom_'+str(central_kpc_one_side)+'_kpc'+\
            '.png'
    if self.verbose: 
      print 'save to'
      print '  ',out_f
    plt.savefig(out_f)


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

    self.plot_range_set_by_camera()
    self.plot_all_quant_zoom(central_kpc_one_side, annotate_clump,
                             clump_list_filename)

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
  def __init__(self, Q_gas, Q_star):

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

    self.debug = Q_gas.debug or Q_star.debug
    self.Q_star_val = Q_star.Q

    self.Q_gas = Q_gas
    self.Q_star = Q_star

    assert self.Q_gas.plane == self.Q_star.plane
    self.plane = self.Q_gas.plane

    assert self.Q_gas.isnap == self.Q_star.isnap
    self.isnap = self.Q_gas.isnap

    self.interpolate_gas_onto_star_grid()


  def interpolate_gas_onto_star_grid(self):
    from scipy import interpolate

    xx = np.linspace(self.Q_gas.coords[self.plane][0].value.min(), self.Q_gas.coords[self.plane][0].value.max(), len(self.Q_gas.coords[self.plane][0]))
    yy = np.linspace(self.Q_gas.coords[self.plane][1].value.min(), self.Q_gas.coords[self.plane][1].value.max(), len(self.Q_gas.coords[self.plane][1]))

    f = interpolate.interp2d(xx,
                             yy,
                             self.Q_gas.Q,
                             kind='cubic')
    Q_gas_resampled = f(xx, yy)
    assert self.Q_star_val.shape == Q_gas_resampled.shape

    if self.debug:
      plt.figure()
      plt.subplot(121)
      plt.imshow(np.log10(self.Q_gas.Q))
      plt.colorbar()
      plt.title('original Qgas')
      plt.subplot(122)
      plt.imshow(np.log10(Q_gas_resampled))
      plt.colorbar()
      plt.title('resmapled Qgas')
      plt.show(block=False)

    self.Q_gas_val = Q_gas_resampled


  def compute_T(self, veldisp_vert, veldisp_r):
    """ see Eqn of Inoue+16"""

    if self.verbose: 
      print "vel disp ratio (sigma_z/sigma_r): ", veldisp_vert / veldisp_r

    res1 = 1. + 0.6 * (veldisp_vert / veldisp_r)**2
    res2 = 0.8 * 0.7 * (veldisp_vert / veldisp_r)
    res = np.where(veldisp_vert < 0.5 * veldisp_r, res1, res2)
    return res

  def compute_T_s(self):
    self.T_s = self.compute_T(self.Q_star.veldisp_vertical_plane.value, self.Q_star.sigma_r)


  def compute_T_g(self):
    self.T_g = self.compute_T(self.Q_gas.veldisp_vertical_plane.value, self.Q_gas.sigma_r)


  def calc_Q_eff(self):

    # The effect of thickness is to increase the stability parameter of each
    # component by a factor T, which depends on the ratio of vertical to
    # radial velocity dispersion.

    w = 2. * self.Q_star.sigma_r * self.Q_gas.sigma_r / (self.Q_star.sigma_r**2 + self.Q_gas.sigma_r**2)

    # 2D array
    res1 = w / (self.Q_star_val * self.T_s) + 1 / (self.Q_gas_val * self.T_g)
    res2 = 1 / (self.Q_star_val * self.T_s) + w / (self.Q_gas_val * self.T_g)

    Q_twoComp_inv = np.where(self.T_s * self.Q_star_val >= self.T_g * self.Q_gas_val, res1, res2)

    self.Q_twoComp = 1. / Q_twoComp_inv
    if self.megaverbose: 
      print(np.isnan(self.Q_twoComp) == True).any()

    return self.Q_twoComp


  def plot_Q_eff(self):
    plt.figure()
    im = plt.imshow(np.log10(self.Q_twoComp),
                    extent=(self.Q_gas._xmin, self.Q_gas._xmax, self.Q_gas._ymin, self.Q_gas._ymax),
                    cmap=cmap_div,
                    origin='lower',
                    vmin=-1, vmax=1)     # clip at -1 < log10(Q) < 1

    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q_{\rm eff}}$", fontsize=16)
    plt.title(r'$Q_{\rm eff}$')
    plt.xlabel('kpc', fontsize=16)
    plt.ylabel('kpc', fontsize=16)

    plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_toomreEff_proj_' + self.plane + '.png'
    if self.verbose: 
      print 'save to'
      print '  ',out_f
    plt.savefig(out_f)


  def plot_Q_eff_zoom(self, central_kpc_one_side=None, annotate_clump=False, clump_list_filename=None):

    if not central_kpc_one_side:
      central_kpc_one_side = 1.5

    if annotate_clump:
      assert clump_list_filename is not None

    xspacing = (self.Q_gas._xmax - self.Q_gas._xmin)/len(self.Q_twoComp)
    xruler = np.arange(self.Q_gas._xmin, self.Q_gas._xmax, xspacing)
    rightBound = np.argmin(abs(xruler - central_kpc_one_side))
    leftBound = np.argmin(abs(xruler + central_kpc_one_side))

    yspacing = (self.Q_gas._ymax - self.Q_gas._ymin)/len(self.Q_twoComp)
    yruler = np.arange(self.Q_gas._ymin, self.Q_gas._ymax, yspacing)
    topBound = np.argmin(abs(yruler - central_kpc_one_side))
    bottomBound = np.argmin(abs(yruler + central_kpc_one_side))

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
        plt.plot(cly, clz, 'x', markersize=15, color='darkturquoise')
      elif self.plane == '1':
        plt.plot(clx, clz, 'x', markersize=15, color='antiquewhite')
      elif self.plane == '2':
        plt.plot(clx, cly, 'x', markersize=15, color='antiquewhite')

    cbar = plt.colorbar(im, extend='both',    # arrows in both direction
                         ticks=[-1, 0, 1]
                        )
    cbar.ax.set_yticklabels([r'$<-1$', r'$0$', r'$>1$'])
    cbar.set_label(r"$\log{Q_{\rm eff}}$", fontsize=16)
    plt.title(r'$Q_{\rm eff}$')
    plt.xlabel('kpc', fontsize=16)
    plt.ylabel('kpc', fontsize=16)

    plt.show(block=False)
    out_f = 'ss' + str(self.isnap) + '_toomreEff_proj_' + self.plane + 'zoomed.png'
    if self.verbose: 
      print 'save to'
      print '  ',out_f
    plt.savefig(out_f)


  def run(self, central_kpc_one_side, annotate_clump, clump_list_filename):

    self.compute_T_g()
    self.compute_T_s()
    self.calc_Q_eff()

    self.plot_Q_eff_zoom(central_kpc_one_side, annotate_clump, clump_list_filename)


if __name__ == '__main__':

  plane     = '0'
  isnap     = 28
  annotate  = False

  testfile  = 'ss'+str(isnap)+'_h2density_clumppos_ncut_0.32_Ncellmin_10.txt'

  Q_gas_obj = ToomreAnalyze(isnap=isnap, wg_var='density',
                      field_type='gas', plane=plane,
                      bin_size_in_log10=0.35, debug=False)
  Q_gas_val = Q_gas_obj.run(radial_nbins=100, central_kpc_one_side=1.5,
                             annotate_clump=annotate,
                             clump_list_filename=testfile)
  Q_gas_obj.plot_all_quant_zoom(1.0, annotate_clump=annotate)

  # # something about calc_kappa doesn't work for stellar component...
  # Q_star_obj = ToomreAnalyze(isnap=isnap, wg_var='mass',
  #                       field_type='star', plane=plane,
  #                       bin_size_in_log10=0.35, debug=True)
  # Q_star_val = Q_star_obj.run(radial_nbins=100, central_kpc_one_side=1.5,
                             #   annotate_clump=True,
                             # clump_list_filename=testfile)

  # Q_tot_obj = ToomreAnalyze_2comp(Q_gas_obj, Q_star_obj)
  # Q_tot_val = Q_tot_obj.run(annotate_clump=True,
#                             clump_list_filename=testfile)


