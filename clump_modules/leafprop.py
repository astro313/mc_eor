'''

Hard coding young stars defined as 10.0 Myr, old stars defined as 100.0 Myr.

Cloud Object to store cloud physical properties.


last mod: 24 July 2018

NOTE
----
- surface area is claculated simply using pi*R^2, perhaps we need to project along different axes given that the clumps/leafs identified are not spherical..?


'''

# print (__doc__)

import numpy as np
import cPickle as pickle
import astropy.constants as C
import matplotlib.pyplot as plt


class Cloud(object):

    def __init__(self, dx, leaf_0_dict, cl_indx):

        self.cl_indx = cl_indx

        self.pc2cm = 3.086E+18
        self.km2cm = 1.E+05
        self.m2cm = 100.0

        self.kg2g = 1.E+03
        self.Msun2g = 1.989E+33
        self.g2Msun = 1./self.Msun2g
        self.kg2Msun = 1. / 1.989E30

        self.s2Myr = 3.17098e-8 * 1.e-6

        self.G_cgs = C.G.value * (self.m2cm)**3 / self.kg2g
        self.J2erg   = 1.E+07
        self.k_B_erg = C.k_B.value * self.J2erg

        self.dx = float(dx)            # pc
        # for the units, see fetch_ga_fields.py
        self.density = leaf_0_dict['density']      # 1/cc
        self.avg_density = np.mean(self.density)   # 1/cc
        self.avg_density_Msun = self.avg_density * C.m_p * self.kg2Msun
        self.mean_density_mass_avg = np.sum(self.density *self.density )/np.sum(self.density)   # 1/cc
        #print self.avg_density_Msun
        self.H2density = leaf_0_dict['H2density']  # 1/cc
        self.Pressure = leaf_0_dict['Pressure']    # K cm-3
        self.P_nt = leaf_0_dict['P_nt']            # K cm-3
        self.metallicity = leaf_0_dict['metallicity']   # old solar units
        self.velx = leaf_0_dict['velx_gas']            # km/s
        self.vely = leaf_0_dict['vely_gas']
        self.velz = leaf_0_dict['velz_gas']
        # plt.hist(self.velx)
        # plt.show()
        # import pdb; pdb.set_trace()

        self.mstar = leaf_0_dict['mass']           # Msun
        self.epoch = leaf_0_dict['epochMyr']       # Myr (age of universe when SF)
        self.mstar_young = leaf_0_dict['young10']       # Myr (age of universe when SF)
        self.mstar_old = leaf_0_dict['old100']       # Myr (age of universe when SF)
        self.velx_star = leaf_0_dict['velx']
        self.vely_star = leaf_0_dict['vely']
        self.velz_star = leaf_0_dict['velz']

        self.vol_pc3()
        self.vol_cc()
        self.mass_Msun()
        self.spherical_radius()

        self.temp()
        self.JeansL_pc()

        self.tff()
        self.massSD()
        self.Mjeans()
        self.Mach_vec()

        self.mass_star_Msun()
        self.young_SFR()
        self.old_SFR()

        self.mean_sigma_NT_mass_avg()
        self.mean_veldisp_mass_avg()
        self.mean_veldisp_mass_avg_stars()
        self.alpha_vir_summed()


    def mean_veldisp_mass_avg_stars(self):

        # mass weighted stellar velocity dispersion
        # cm/s

        mean  = []
        for vec in [self.velx_star,self.vely_star,self.velz_star]:
          xx    = np.sum(self.mstar*vec)/np.sum(self.mstar)
          mean.append(xx)

        std  = []
        for ii,vec in enumerate([self.velx_star,self.vely_star,self.velz_star]):
          xx   = np.sum(self.mstar*(vec[:]-mean[ii])**2)/np.sum(self.mstar)
          std.append(1.e+5*np.sqrt(xx))

        self.mean_veldisp_mass_avg_stars = np.array(std)


    def mean_veldisp_mass_avg(self):

        # mass weighted bulk gas velocity dispersion
        # cm/s

        mean  = []
        for vec in [self.velx,self.vely,self.velz]:
          xx    = np.sum(self.density*vec)/np.sum(self.density)
          mean.append(xx)

        std  = []
        for ii,vec in enumerate([self.velx,self.vely,self.velz]):
          xx   = np.sum(self.density*(vec[:]-mean[ii])**2)/np.sum(self.density)
          std.append(1.e+5*np.sqrt(xx))

        self.mean_veldisp_mass_avg = np.array(std)

    def mean_sigma_NT_mass_avg(self):
        x     =  (self.P_nt * self.k_B_erg / C.m_p.cgs.value) / self.density
        x     =  np.sqrt(x)
        self.mean_sigma_NT_mass_avg =  np.sum(self.density * x) / np.sum(self.density)

    def alpha_vir_summed(self, include_bulk=True, include_starVel=True):

      x = (self.Pressure * self.k_B_erg / C.m_p.cgs.value) / self.density
      self.cs_avg = np.sqrt(np.sum(self.density * x)/np.sum(self.density))

      v_NT = self.mean_sigma_NT_mass_avg/1.e+5
      v_cs = self.cs_avg/1.e5

      if include_bulk:
        v_disp = 0.0
        for i,x in enumerate(['x','y','z']):
          v_disp = v_disp + (self.mean_veldisp_mass_avg[i]/1.e+5)**2
        v_disp = (1.0/3.0)*np.sqrt(v_disp)

        sigma_2_tot = v_NT**2 + v_cs**2 + v_disp**2
      else:
        sigma_2_tot = v_NT**2 + v_cs**2

      if include_starVel:
        v_disp_stars = 0.0
        for i,x in enumerate(['x','y','z']):
          v_disp_stars = v_disp_stars + (self.mean_veldisp_mass_avg_stars[i]/1.e+5)**2
        v_disp_stars = (1.0/3.0)*np.sqrt(v_disp_stars)

        sigma_2_tot  = (self.mass_Msun*sigma_2_tot + self.mstar_Msun_tot*v_disp_stars**2)/(self.mstar_Msun_tot + self.mass_Msun)

      from astropy import units,constants
      conv = (units.km/units.s)**2 * constants.pc /( constants.G * constants.M_sun)
      conv = conv.to(1).value

      mass = self.mstar_Msun_tot + self.mass_Msun

      self.alpha_vir_summed  = 5.0 * sigma_2_tot * conv * self.R_pc/mass


    # --- star particles ---
    def mass_star_Msun(self):
        self.mstar_Msun_tot = self.mstar.sum()
        self.s2gR = self.mstar_Msun_tot/self.mass_Msun

    def young_SFR(self):
        """ calculate the SFR based on young stars [Msun/yr]
        """
        self.young_SFR_MsunPyr = self.mstar_young.sum() / 10.0e6

    def old_SFR(self):
        """ calculate the SFR based on old stars [Msun/yr]
        """
        self.old_SFR_MsunPyr = self.mstar_old.sum() / 100.0e6

    # --- AMR ----
    def vol_pc3(self):
        self.vol_pc3 = self.dx**3 * len(self.density)

    def vol_cc(self):
        self.vol_cc = (self.dx * self.pc2cm)**3

    def mass_Msun(self):
        # in unit of Msun
        self.mass_Msun = (self.density * self.vol_cc).sum() * \
            C.m_p.value * self.kg2Msun

    def spherical_radius(self):
        """ assume spherical cloud. """
        self.R_pc = (3. * self.vol_pc3 / 4. / np.pi)**(1 / 3.)
        self.R_cm = self.R_pc * self.pc2cm

    def temp(self):   # get average temp??? mass-weighted???
        """ Temp in Kelvin. """
        self.T = self.Pressure / self.density

    def JeansL_pc(self, total_pressure= True):

        temp   = self.T
        if total_pressure:
          temp = temp + self.P_nt/self.density

        self.L_jeans_pc = 3.31 * \
            (self.density / 100.0)**(-0.5) * (temp / 20.0)**0.5

    def tff(self):
        self.tff = (3. * np.pi / 32.0 / self.G_cgs / (self.avg_density * C.m_p.value * self.kg2g))**0.5
        self.tff_Myr = self.tff * self.s2Myr


    def massSD(self):
        """ assuming spherical clouds --> circular in projection. """
        self.massSD = self.mass_Msun / (np.pi * self.R_pc**2)


    def Mjeans(self):
        """

        Cloud's Jeans Mass as computed in Joung & Mac Low 2006

        Mj = rho_avg * lambdaJeans ^3

        lambdaJeans = (pi / G rhoavg)^1/2 sigma_tot
        sigma_tot   = (cs^2 + 1/3 sigma3D^2l ) ^1/2


        """

        self.M_jeans = (self.avg_density * C.m_p.cgs.value) * \
            (np.average(self.L_jeans_pc) * self.pc2cm)**3 * self.g2Msun


    def Mach_vec(self):
        """

        Calculate from mass-weighted NT and thermal pressure.

        """
        mach = np.sqrt(1 +  self.P_nt/self.Pressure)
        self.Mach_vec = np.sqrt((self.density * mach).sum()/self.density.sum())


    def __str__(self):
        print '\n', '=' * 100
        print("Calculated parameters of cloud  {:d}").format(self.cl_indx)
        print("Gas mass x 10^7              = {:.2f} [Msun]").format(self.mass_Msun/1.e7)
        print("Stellar mass     = {:.2f} x 10^7 [Msun]").format(self.mstar_Msun_tot/1.e7)
        print("Volume                   = {:.2f} [pc^3]").format(self.vol_pc3)
        print("Number of Cells          = {:d}").format(len(self.density))
        print("average density          = {:.2f}").format(self.avg_density)
        print("Spherical radius         = {:.2f} [pc]").format(self.R_pc)
        print "Mass Surface Density     = {:.2f} [Msun/pc^2]".format(self.massSD)
        print("Free fall time           = {:.2f} [Myr]").format(self.tff_Myr)

        print "*" * 10 + " mass weighted averages: " + "*" * 10
        print(" density          = {:.2f} cm-3").format(self.mean_density_mass_avg)           # mass weighted density
        print(" Mach             = {:.2f} ").format((np.mean(self.Mach_vec)))                 # as in vallini+18
        print(" v turb           = {:.2f} km/s").format(self.mean_sigma_NT_mass_avg/1.e+5)    # velocity from  non thermal pressure component
        v_disp = 0.0
        for i,x in enumerate(['x','y','z']):
          print(" v disperison {}  = {:.2f} km/s").format(x,self.mean_veldisp_mass_avg[i]/1.e+5) # 1d bulk motion mass weighted std(v)
          v_disp = v_disp + (self.mean_veldisp_mass_avg[i]/1.e+5)**2
        v_disp = (1.0/3.0)*np.sqrt(v_disp)
        v_disp_stars = 0.0
        for i,x in enumerate(['x','y','z']):
          print(" v disp stars {}  = {:.2f} km/s").format(x,self.mean_veldisp_mass_avg_stars[i]/1.e+5) # 1d bulk motion mass weighted std(v)
          v_disp_stars = v_disp_stars + (self.mean_veldisp_mass_avg_stars[i]/1.e+5)**2
        v_disp_stars = (1.0/3.0)*np.sqrt(v_disp_stars)
        print(" bulk dispersion stars  = {:.2f} km/s").format(v_disp_stars)

        print(" cs               = {:.2f} [km/s]").format(self.cs_avg/1.e5)
        print(" alpha            = {:.2f} ").format(self.alpha_vir_summed)

        return '=' * 100



# below is for testing
if __name__ == '__main__':

    import sys
    sys.path.append('../')

    from io_modules.manipulate_fetch_gal_fields import get_units, get_dx
    import pymses

    snapshot_num = 16
    leafdir = 'test_brute/leaf_fields_' + str(snapshot_num) + '/'

    # fname = "0.32_10_fields.p"
    fname = "6.81_10_fields.p"

    leaf_fields = pickle.load(open(leafdir + fname, "rb"))

    for cidx in range(len(leaf_fields)):
        print Cloud(get_dx(snapshot_num), leaf_fields[str(cidx)], cidx)

