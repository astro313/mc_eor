'''

last mod: 16 July 2018

'''


import numpy as np


def import_fetch_gal(isnap=28, folder_ramsesdata='output', tag_h5file="_center_fields0123456-15_resampled.h5", verbose=True):

    import pymses
    from pymses.utils import constants as C
    import h5py
    from fetch_gal_fields import get_units

    nsnap = str(isnap)

    h5_file = "snapshot" + nsnap + tag_h5file
    if verbose:
        print 'reading output from fetch_gal_fields() from'
        print '  ', h5_file

    f = h5py.File(h5_file, "r")
    # careful sometimes i used "density" (e.g., resample.py), see
    # resample_fields.py to make sure
    density = f["rho"].value
    H2 = f["H2"].value
    Pressure = f["P"].value
    P_nt = f["P_nt"].value
    metallicity = f["Z"].value
    velx = f["vel_x"].value
    vely = f["vel_y"].value
    velz = f["vel_z"].value

    ro = pymses.RamsesOutput(folder_ramsesdata, isnap)

    factor_density, unit_dens = get_units(ro=ro)['rho']      # 1/cm^3 (not H/cm^3)
    density *= factor_density
    if(verbose):
        print 'max density'
        print density.max(), unit_dens

    factor_vel, unit_vel = get_units(ro=ro)['vel']
    velx *= factor_vel
    vely *= factor_vel
    velz *= factor_vel
    if(verbose):
        print 'max vel'
        print velx.max(), vely.max(), velz.max(), unit_vel

    factor_P, unit_P = get_units(ro=ro)['P']
    Pressure *= factor_P
    P_nt *= factor_P
    if(verbose):
        print 'max P, P_nt'
        print np.log10(Pressure.max()), np.log10(P_nt.max()), unit_P

    data = dict(density=density, H2=H2,
                P=Pressure,
                P_nt=P_nt,
                Z=metallicity,
                velx=velx,
                vely=vely,
                velz=velz
                )
    if True:
        if verbose:
            print 'Clipping variables'
        data["density"][data["density"] <= 0] = np.min(
            data["density"][data["density"] > 0])
        data["H2"][data["H2"] < 1.e-3] = 1.e-3
    if verbose:
        for var in ['density', 'H2']:
            print '  ', var, np.max(data[var]), np.min(data[var])

    return data


def prepare_unigrid(data, verbose=False):

    import yt

    if(not verbose):
        from yt.funcs import mylog
        mylog.setLevel(40)

    field = ("h2density")

    def _h2density(field, data):
        try:
            return data["density"] * data["H2"]
        except:
            return data[("stream", "density")] * data[("stream", "H2")]

    ds = yt.load_uniform_grid(data, data["density"].shape)
    dd = ds.all_data()
    ds.add_field(("stream", "h2density"), function=_h2density, units="g/cm**3")  # unit is in g/cc only if convert_unit is properly called when loading in data
    assert (dd['H2'] * dd['density']).max() == dd['h2density'].max()

    return ds, dd


def check_hist_h2(data):
    import matplotlib.pyplot as plt
    aa = (data["density"] * data["H2"]).flatten()
    print np.max(aa), np.min(aa), np.min(aa[aa > 0])
    aa[aa <= 0] = np.min(aa[aa > 0])
    aa = np.log10(aa)

    print np.max(aa), np.min(aa), np.min(aa[aa > 0])

    plt.close('all')
    plt.figure()
    plt.hist(aa, bins=100)
    for ele in th_list:
        x = np.log10(ele)
        plt.plot([x, x], [1, 1.e+7], ls='--', color='k')
    plt.yscale('log')
    plt.savefig('hist_test.png')
