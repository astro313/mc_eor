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
        self.spherical_radius_pc()
        self.spherical_radius_cm()
        self.temp()
        self.JeansL_pc()
        self.veldisp()
        self.sigmaSq_3()
        self.sigmaSq()
        self.sigmaSq_2()
        self.sigmaSq_NT()
        self.tot_veldisp()
        self.tff()
        self.tff_Myr()
        self.massSD()
        self.turb_KE()
        self.Mjeans()
        self.alphaVir()
        self.Mach_vec()
        self.SFR_per_tff()

        self.mass_star_Msun()
        self.young_SFR()
        self.old_SFR()
        self.alphaVir_total()

        self.mean_sigma_NT_mass_avg()
        self.mean_veldisp_mass_avg()
        self.alpha_vir_summed()

    def tot_veldisp(self):
        """ non-thermal, turbulent + bulk + average mass-weighted sound speed"""
        x = (self.Pressure * self.k_B_erg / C.m_p.cgs.value) / self.density
        self.cs_avg = np.sqrt(np.sum(self.density * x)/np.sum(self.density))

        self.sigmaSq_bulk_cs = self.sigmaSq + self.cs_avg**2
        self.sigmaSq_tot = self.sigmaSq_bulk_cs + self.sigmaSq_NT

        #print "sigma^2 total before and including pressure: "
        #print np.sqrt(self.sigmaSq_bulk_cs)/1.e5, np.sqrt(self.sigmaSq_tot)/1.e5


    def alphaVir(self):
        """

        Cloud's alpha virial, Bertoldi & McKee 1992.

                   5 * sigma ^ 2 * R
          alpha =  ------------------
                        G * M

        """

        # bulk velo and cs only, excluding NT pressure
        self.alpha = 5. * (self.sigmaSq_bulk_cs) * \
            self.R_cm / (self.G_cgs * self.mass_Msun * self.Msun2g)

        # sigma from NT and cs only, excluding bulk
        self.alpha_NT = 5. * (self.sigmaSq_NT + self.cs_avg**2) * \
            self.R_cm / (self.G_cgs * self.mass_Msun * self.Msun2g)

        # sigma from all terms
        self.alpha_bulk_NT = 5. * (self.sigmaSq_tot) * \
            self.R_cm / (self.G_cgs * self.mass_Msun * self.Msun2g)

    def alpha_vir_summed(self):

      v_NT = self.mean_sigma_NT_mass_avg/1.e+5
      v_cs = self.cs_avg/1.e5

      if True:
        v_disp = 0.0
        for i,x in enumerate(['x','y','z']):
          v_disp = v_disp + (self.mean_veldisp_mass_avg[i]/1.e+5)**2
        v_disp = (1.0/3.0)*np.sqrt(v_disp)

        sigma_2_tot = v_NT**2 + v_cs**2 + v_disp**2
      else:
        sigma_2_tot = v_NT**2 + v_cs**2


      from astropy import units,constants
      conv = (units.km/units.s)**2 * constants.pc /( constants.G * constants.M_sun)
      conv = conv.to(1).value

      mass = self.mstar_Msun_tot + self.mass_Msun

      self.alpha_vir_summed  = 5.0 * sigma_2_tot * conv * self.R_pc/mass


    def alphaVir_total(self):
        """

        Cloud's alpha virial, Bertoldi & McKee 1992.

                   5 * sigma ^ 2 * R
          alpha =  ------------------
                        G * M

        """

        # mass-weighted stellar velo disp
        _velx_star = (self.velx_star - np.mean(self.velx_star)) * self.km2cm
        _vely_star = (self.vely_star - np.mean(self.vely_star)) * self.km2cm
        _velz_star = (self.velz_star - np.mean(self.velz_star)) * self.km2cm

        self.sigmaSq_star = 1. / 3 * \
            (self.mstar * ((_velx_star)**2 + (_vely_star) **
                             2 + (_velz_star)**2)).sum() / self.mstar.sum()

        #print "stellar sigma [km/s]: ", np.sqrt(self.sigmaSq_star) / 1.e5    # km/s
        del _velx_star, _vely_star, _velz_star

        self.alpha_total = 5. * (self.mass_Msun * self.Msun2g * self.sigmaSq_bulk_cs + self.mstar_Msun_tot *self.Msun2g * self.sigmaSq_star) / self.G_cgs / ((self.mass_Msun * self.Msun2g)**2 / self.R_cm + (self.mstar_Msun_tot * self.Msun2g)**2 / self.R_cm)

        self.alpha_NT_total = 5. * (self.mass_Msun * self.Msun2g * (self.sigmaSq_NT + self.cs_avg**2) + self.mstar_Msun_tot *self.Msun2g * self.sigmaSq_star) / self.G_cgs / ((self.mass_Msun * self.Msun2g)**2 / self.R_cm + (self.mstar_Msun_tot * self.Msun2g)**2 / self.R_cm)

        self.alpha_all_total = 5. * (self.mass_Msun * self.Msun2g * self.sigmaSq_tot + self.mstar_Msun_tot *self.Msun2g * self.sigmaSq_star) / self.G_cgs / ((self.mass_Msun * self.Msun2g)**2 / self.R_cm + (self.mstar_Msun_tot * self.Msun2g)**2 / self.R_cm)

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

    def spherical_radius_pc(self):
        """ assume spherical cloud. """
        self.R_pc = (3. * self.vol_pc3 / 4. / np.pi)**(1 / 3.)

    def spherical_radius_cm(self):
        """ assume spherical cloud. """
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

    def veldisp(self):
        self.xdisp = np.std(self.velx) * self.km2cm
        self.ydisp = np.std(self.vely) * self.km2cm
        self.zdisp = np.std(self.velz) * self.km2cm

        # plt.hist(self.velx, bins=30)
        # plt.show()

    def sigmaSq_3(self):
        sss = 1. / 3 * (self.density * (self.xdisp**2 + self.ydisp **
                                        2 + self.zdisp**2)).sum() / self.density.sum()
        #print "sigma calculated slightly different [km/s]: ", np.sqrt(sss) / 1.e5

    def sigmaSq(self):
        """ mass-weighted 1-D sigmaSq """

        _velx = (self.velx - np.mean(self.velx)) * self.km2cm
        _vely = (self.vely - np.mean(self.vely)) * self.km2cm
        _velz = (self.velz - np.mean(self.velz)) * self.km2cm

        self.sigmaSq = 1. / 3 * \
            (self.density * ((_velx)**2 + (_vely) **
                             2 + (_velz)**2)).sum() / self.density.sum()

        #print "sigma from vmean = np.mean() [km/s]: ", np.sqrt(self.sigmaSq) / 1.e5    # km/s
        del _velx, _vely, _velz

    def sigmaSq_2(self):
        """ mass-weighted 1-D sigmaSq """

        _x_mean = (self.density * (self.velx)).sum() / self.density.sum()
        _y_mean = (self.density * (self.vely)).sum() / self.density.sum()
        _z_mean = (self.density * (self.velz)).sum() / self.density.sum()

        _velx = (self.velx - _x_mean) * self.km2cm
        _vely = (self.vely - _y_mean) * self.km2cm
        _velz = (self.velz - _z_mean) * self.km2cm

        sigmaSq = 1. / 3 * (self.density * ((_velx)**2 +
                                            (_vely)**2 +
                                            (_velz)**2)).sum() / self.density.sum()
        # km/s
        #print "sigma calculated with mass-weighted vmean [km/s]: ", np.sqrt(sigmaSq) / 1.e5
        del _x_mean, _y_mean, _z_mean, _velx, _vely, _velz


    def sigmaSq_NT(self):
        x = (self.P_nt * self.k_B_erg / C.m_p.cgs.value) / self.density
        self.sigmaSq_NT = np.sum(self.density * x) / np.sum(self.density)

    def mean_veldisp_mass_avg(self):

        # mass weighted velocity dispersion
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

    def tff(self):
        self.tff = (3. * np.pi / 32.0 / self.G_cgs / (self.avg_density * C.m_p.value * self.kg2g))**0.5

    def tff_Myr(self):
        self.tff_Myr = self.tff * self.s2Myr

    def massSD(self):
        """ assuming spherical clouds --> circular in projection. """
        self.massSD = self.mass_Msun / (np.pi * self.R_pc**2)

    def turb_KE(self):
        self.turb_KE = 1 / 2. * self.mass_Msun * self.sigmaSq

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
        #print "Mass-weighted Mach from pressure: "
        mach = np.sqrt(1 +  self.P_nt/self.Pressure)
        self.Mach_vec = np.sqrt((self.density * mach).sum()/self.density.sum())
        #print self.Mach_vec


    def SFR_per_tff(self):
        """ Find the star formation  per free fall time of the cloud. As in  Krumholz, Matzner & McKee 2006. Eqn 40 and 41.

                                     -0.68        -0.32
        sfr_ff = 0.073 * alpha_vir ^       Mach ^
        SFR = Mdot = sfr_ff * M_cl / t_ff

        """
        self.sfr_ff = 0.073 * self.alpha**(-0.68) * np.mean(self.Mach_vec)**(-0.32)
        self.SFR = self.sfr_ff * self.mass_Msun / self.tff_Myr

        # Star formation rate as in Joung & Mac Low 2006.
        if (self.mass_Msun / self.M_jeans > 1.0):
            self.SFR_JML = 0.3 * self.mass_Msun / self.tff_Myr
        else:
            self.SFR_JML = 0.0

    def __str__(self):
        print '\n', '=' * 100
        print("Calculated parameters of cloud  {:d}").format(self.cl_indx)
        print("Mass x 10^7              = {:.2f} [Msun]").format(self.mass_Msun/1.e7)
        print("Volume                   = {:.2f} [pc^3]").format(self.vol_pc3)
        print("Number of Cells          = {:d}").format(len(self.density))
        print("average density          = {:.2f}").format(self.avg_density)
        print("Spherical radius         = {:.2f} [pc]").format(self.R_pc)
        print "Mass Surface Density     = {:.2f} [Msun/pc^2]".format(self.massSD)
        print("Free fall time           = {:.2f} [Myr]").format(self.tff_Myr)
        print("Velocity disp            = {:.2f}, {:.2f}, {:.2f} [km/s]").format(self.xdisp/1.e5, self.ydisp/1.e5, self.zdisp/1.e5)
        print("Velocity disp 3D         = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq)/1.e5)

        """

        print "*" * 10 + "   Jeans Mass Calculation:   " + "*" * 10
        print("cs avg       = {:.2f} [km/s]").format(self.cs_avg/1.e5)
        print("sigma from velo   = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq)/1.e5)
        print("sigma from NT Pressure = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq_NT)/1.e5)
        print("sigma tot before adding NT Pressure = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq_bulk_cs)/1.e5)
        print("sigma tot    = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq_tot)/1.e5)
        print("Jeans length = {:.2f} [pc]").format(np.mean(self.L_jeans_pc))
        print("Jeans Mass   = {:.2f} [Msun]").format(self.M_jeans)

        print "*" * 10 + "   Star Formation rate calculation:   " + "*" * 10
        print("alpha virial from bulk velo + cs only = {:.2f}").format(self.alpha)
        print("alpha virial from bulk velo incl. stellar  = {:.2f}").format(self.alpha_total)
        print("alpha virial from NT velo + cs only = {:.2f}").format(self.alpha_NT)
        print("alpha virial from NT velo + incl. stellar  = {:.2f}").format(self.alpha_NT_total)
        print("alpha virial from all sigma, exclude star = {:.2f}").format(self.alpha_bulk_NT)
        print("alpha virial from all sigma, incl. star = {:.2f}").format(self.alpha_all_total)
        print("Mach number from NT pressure = {:.2f} ").format((np.mean(self.Mach_vec)))
        print("sfr_tff       = {:.2f} ").format(self.sfr_ff)
        print("SFR          = {:.2f} [Msun/yr] ").format(self.SFR / 1.e6)
        print("SFR simple 30% SFE   = {:.2f} [Msun/yr]").format(self.SFR_JML/1.e6)
        print "*" * 10 + " star particle properties: " + "*" * 10
        print ("Stellar mass: {:.2f} x10^7 Msun").format(self.mstar_Msun_tot/1.e7)
        print ("stellar to gas mass ratio: {:.4f}").format(self.s2gR)
        print ("SFR based on young stars in structure: {:.4f} [Msun/yr] ").format(self.young_SFR_MsunPyr)
        print ("SFR based on old stars in structure: {:.4f} [Msun/yr] ").format(self.old_SFR_MsunPyr)

        """

        print "*" * 10 + " mass weighted averages: " + "*" * 10
        print(" gas     mass     = {:.2f} x 10^7 [Msun]").format(self.mass_Msun/1.e7)         # gas mass
        print(" stellar mass     = {:.2f} x 10^7 [Msun]").format(self.mstar_Msun_tot/1.e7)    # stellar masss
        print(" density          = {:.2f} cm-3").format(self.mean_density_mass_avg)           # mass weighted density
        print(" Mach             = {:.2f} ").format((np.mean(self.Mach_vec)))                 # as in vallini+18
        print(" v turb           = {:.2f} km/s").format(self.mean_sigma_NT_mass_avg/1.e+5)    # velocity from  non thermal pressure component
        v_disp = 0.0
        for i,x in enumerate(['x','y','z']):
          print(" v disperison {}  = {:.2f} km/s").format(x,self.mean_veldisp_mass_avg[i]/1.e+5) # 1d bulk motion mass weighted std(v)
          v_disp = v_disp + (self.mean_veldisp_mass_avg[i]/1.e+5)**2
        v_disp = (1.0/3.0)*np.sqrt(v_disp)
        print(" bulk dispersion  = {:.2f} km/s").format(v_disp)                                   # 3d bulk motion mass weighted std
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

    fname = "0.32_10_fields.p"
    #fname = "6.81_10_fields.p"

    leaf_fields = pickle.load(open(leafdir + fname, "rb"))

    for cidx in range(len(leaf_fields)):
        print Cloud(get_dx(snapshot_num), leaf_fields[str(cidx)], cidx)

