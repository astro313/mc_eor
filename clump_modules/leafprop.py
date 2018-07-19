'''

Cloud Object to store cloud physical properties.

last mod: 18 July 2018

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

        self.dx = float(dx)            # pc
        # for the units, see fetch_ga_fields.py
        self.density = leaf_0_dict['density']      # 1/cc
        self.avg_density = np.mean(self.density)   # 1/cc
        self.avg_density_Msun = self.avg_density * C.m_p * self.kg2Msun
        print self.avg_density_Msun
        self.H2density = leaf_0_dict['H2density']  # 1/cc
        self.Pressure = leaf_0_dict['Pressure']    # K cm-3
        self.P_nt = leaf_0_dict['P_nt']            # K cm-3
        self.metallicity = leaf_0_dict['metallicity']   # old solar units
        self.velx = leaf_0_dict['velx']            # km/s
        self.vely = leaf_0_dict['vely']
        self.velz = leaf_0_dict['velz']
        # plt.hist(self.velx)
        # plt.show()
        # import pdb; pdb.set_trace()
        self.vol_pc3()
        self.vol_cc()
        self.mass_Msun()
        self.spherical_radius_pc()
        self.spherical_radius_cm()
        self.temp()
        self.JeansL_pc()
        self.JeansL_cm()
        self.veldisp()
        self.sigmaSq_3()
        self.sigmaSq()
        self.sigmaSq_2()
        self.tot_veldisp()
        self.tff()
        self.tff_Myr()
        self.massSD()
        self.turb_KE()
        self.Mjeans()
        self.alphaVir()
        self.Mach()
        self.Mach_vec()
        self.SFR_per_tff()

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

    def JeansL_pc(self):
        self.L_jeans_pc = 3.31 * \
            (self.density / 100.0)**(-0.5) * (self.T / 20.0)**0.5

    def JeansL_cm(self):
        J2erg = 1.E+07
        self.k_B_erg = C.k_B.value * J2erg

        mH = C.m_p.value * self.kg2g
        mu = 1.3017 * mH
        self.L_jeans_cm = (15.0 * self.T * self.k_B_erg / 4.0 /
                           np.pi / self.G_cgs / mu**2 / self.density)**0.5

        # # just to make sure the equations we use in JeansL_pc and JeansL_cm are equivalent, which they!
        # plt.figure()
        # plt.hist(self.L_jeans_pc * self.pc2cm)
        # plt.title("110")
        # plt.show(block=False)

        # plt.figure()
        # plt.hist(self.L_jeans_cm)
        # plt.title("116")
        # plt.show(block=False)

    def veldisp(self):
        self.xdisp = np.std(self.velx) * self.km2cm
        self.ydisp = np.std(self.vely) * self.km2cm
        self.zdisp = np.std(self.velz) * self.km2cm

        # plt.hist(self.velx, bins=30)
        # plt.show()

    def sigmaSq_3(self):
        sss = 1. / 3 * (self.density * (self.xdisp**2 + self.ydisp **
                                        2 + self.zdisp**2)).sum() / self.density.sum()
        print "sigma calculated slightly different [km/s]: ", np.sqrt(sss) / 1.e5

    def sigmaSq(self):
        """ mass-weighted 1-D sigmaSq """

        _velx = (self.velx - np.mean(self.velx)) * self.km2cm
        _vely = (self.vely - np.mean(self.vely)) * self.km2cm
        _velz = (self.velz - np.mean(self.velz)) * self.km2cm

        self.sigmaSq = 1. / 3 * \
            (self.density * ((_velx)**2 + (_vely) **
                             2 + (_velz)**2)).sum() / self.density.sum()

        print "sigma [km/s]: ", np.sqrt(self.sigmaSq) / 1.e5    # km/s
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
        print "sigma calculated slightly different [km/s]: ", np.sqrt(sigmaSq) / 1.e5
        del _x_mean, _y_mean, _z_mean, _velx, _vely, _velz

    def tot_veldisp(self):
        """ non-thermal, turbulent + average mass-weighted sound speed"""

        self.cs_avg = np.sqrt(np.mean(self.Pressure * self.k_B_erg / C.m_p.cgs.value) / self.avg_density)
        self.sigmaSq_tot = self.sigmaSq + self.cs_avg**2

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

    def alphaVir(self):
        """

        Cloud's alpha virial, Bertoldi & McKee 1992.

                   5 * sigma ^ 2 * R
          alpha =  ------------------
                        G * M

        """

        self.alpha = 5. * (self.sigmaSq_tot) * \
            self.R_cm / (self.G_cgs * self.mass_Msun * self.Msun2g)

    def Mach(self):
        self.Mach = np.sqrt(self.sigmaSq) / self.cs_avg

    def Mach_vec(self):
        self.Mach_vec = np.sqrt((self.Pressure + self.P_nt)/self.Pressure)
        # weight by density?

        print np.log10(self.Mach), np.log10(np.mean(self.Mach_vec))

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

        print "*" * 10 + "   Jeans Mass Calculation:   " + "*" * 10
        print("cs avg       = {:.2f} [km/s]").format(self.cs_avg/1.e5)
        print("sigma turb   = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq)/1.e5)
        print("sigma tot    = {:.2f} [km/s]").format(np.sqrt(self.sigmaSq_tot)/1.e5)
        print("Jeans length = {:.2f} [pc]").format(np.mean(self.L_jeans_pc))
        print("Jeans Mass   = {:.2f} [Msun]").format(self.M_jeans)

        print "*" * 10 + "   Star Formation rate calculation:   " + "*" * 10
        print("alpha virial  = {:.2f}").format(self.alpha)
        print("Mach number   = {:.2f} ").format(self.Mach)
        print("sfr_tff       = {:.2f} ").format(self.sfr_ff)
        print("SFR          = {:.2f} [Msun/yr] ").format(self.SFR / 1.e6)
        print ""
        print("SFR simple 30% SFE   = {:.2f} [Msun/yr]").format(self.SFR_JML/1.e6)
        return '=' * 100


# below is for testing
if __name__ == '__main__':

    # to load back in fields of each leaf
    snapshot_num = 28
    leafdir = 'leaf_fields_' + str(snapshot_num) + '/'
    # '{0:.2f}'.format(args.ncut) + '_' + str(args.step) + '_' + str(args.Nmin) from snapshot28_leafprop.py
    fname = "0.03_5_3_fields.p"

    leaf_fields = pickle.load(open(leafdir + fname, "rb"))
    outClPickle = fname.replace('.p', '_cloudprop.p')

    # hard-code to get dx for now...
    # saved in fetch_gal_fields.py
    from io_modules.manipulate_fetch_gal_fields import get_units, get_dx
    import pymses

    # call cloud 0 for testing code..
    C0 = Cloud(get_dx(snapshot_num), leaf_fields['0'], int('0'))
    print C0

