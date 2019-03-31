'''

last mod: 16 July 2018

'''

import numpy as np
import pymses


def get_dx(ssnum):
    # hard-code to get dx for now...
    # saved in fetch_gal_fields.py
    ro = pymses.RamsesOutput("output", ssnum)
    ds = np.load('snapshot' + str(ssnum) + '_center_fields0123456-15.npz')
    dx_vector = ds['dx_vector']

    originalLevels = np.log2(1. / np.unique(dx_vector))
    highestRes = 2.**(originalLevels.max() * -1)
    # after finely sample unigrid
    dx_pc = highestRes * get_units(ro=ro)['dx'][0]
    return dx_pc


def get_units(ro=None):
    from pymses.utils import constants as C
    assert ro is not None
    # conversion dictionary
    dict_unit = {}
    dict_unit['dx'] = [ro.info['unit_length'].express(C.pc), 'pc']
    dict_unit['rho'] = [
        (ro.info['unit_density'] / C.mH).express(1 / C.cm**3), 'cm-3']
    dict_unit['P'] = [ro.info['unit_pressure'].express(
        C.erg / C.cm**3) / C.kB.express(C.erg / C.K), 'K cm-3']
    dict_unit['P_nt'] = dict_unit['P']
    dict_unit['H2'] = [1, '']
    dict_unit['vel'] = [ro.info['unit_velocity'].express(C.km / C.s), 'km/s']
    dict_unit['mass'] = [ro.info['unit_mass'].express(C.Msun), 'Msun']
    dict_unit['epoch'] = [ro.info['unit_time'].express(C.Gyr), 'Gyr']

    return dict_unit


def import_fetch_gal(isnap=28, folder_ramsesdata='output', tag_h5file="_center_fields0123456-15_resampled.h5", verbose=True, convert=True, clipping = 'min'):

    import pymses
    from pymses.utils import constants as C
    import h5py

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

    if convert:
        factor_density, unit_dens = get_units(
            ro=ro)['rho']      # 1/cm^3 (not H/cm^3)
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
    if clipping is not None:

      if clipping == 'min':
        clip_n  = np.min(data["density"][data["density"] > 0])
        clip_H2 = 1.e-3
      else:
        clip_n,clip_H2  = clipping[0],clipping[1]

        if verbose:
            print 'Clipping variables'
        data["density"][data["density"] <= clip_n ] = clip_n
        data["H2"     ][data["H2"]      <= clip_H2] = clip_H2
    if verbose:
        for var in ['density', 'H2']:
            print '  ', var, np.max(data[var]), np.min(data[var])

    return data


def import_fetch_stars(isnap=28, folder_ramsesdata='output', tag_h5file="_center_stars_resampled.h5",
                       verbose=True, convert=True, clipping = 'min'):
    """

    In the input .h5 file, epochMyr is converted into Myr (Universe Age).

    Returns
    -------
    dataDict: dict
        containing star fields (mass, epochMyr, young10, old100) in converted units

    """

    import pymses
    from pymses.utils import constants as C
    import h5py

    nsnap = str(isnap)

    h5_file = "snapshot" + nsnap + tag_h5file
    ro = pymses.RamsesOutput(folder_ramsesdata, isnap)

    if verbose:
        print 'reading output from call_resampled_fields.py from'
        print '  ', h5_file

    f = h5py.File(h5_file, "r")

    for fff in f.iterkeys():
        print fff
        exec(str(fff) + '= f[fff].value')

    # print mass.shape, epochMyr.shape, young10.shape, old100.shape, velx.shape

    if verbose:
        print 'max and min of Universe age when stars are formed [Gyr] '
        print epochMyr.max() / 1.e3, epochMyr[epochMyr > 0].min() / 1.e3

    if convert:
        factor_mass, unit_mass = get_units(ro=ro)['mass']
        factor_vel, unit_vel = get_units(ro=ro)['vel']
        mass *= factor_mass
        velx *= factor_vel
        vely *= factor_vel
        velz *= factor_vel

        if(verbose):
            print 'max vel'
            print velx.max(), vely.max(), velz.max(), unit_vel

        for fff in f.iterkeys():
            if ('young' in fff) or ('old' in fff):
                exec(str(fff) + '= factor_mass * f[fff].value')

        if(verbose):
            print 'max mass in log10'
            print np.log10(mass.max()), unit_mass
            print ''
            print 'sum of all masses in log: Msun'
            print np.log10(mass.sum())

    dataDict = {}
    for fff in f.iterkeys():
        dataDict[fff] = eval(fff)

    if clipping is not None:
      if clipping == 'min':
          clip_val = np.min(dataDict["mass"][dataDict["mass"] > 0])
      else:
          clip_val = clipping
      if verbose:
        print 'Clipping variables at',clip_val

      dataDict["mass"][dataDict["mass"] <= clip_val] = clip_val

      if verbose:
        for var in ['mass']:
            print ' Mass range in log10, after clipping: ', var, np.log10(np.max(dataDict[var])), np.log10(np.min(dataDict[var]))

    return dataDict


def prepare_unigrid(data, regionsize_kpc=7., verbose=False, add_unit=False, debug=False):

    import yt

    if(not verbose):
        from yt.funcs import mylog
        mylog.setLevel(40)

    if(add_unit):
        # see http://yt-project.org/doc/examining/generic_array_data.html
        # for reference

        from astropy import constants as cc
        mp = cc.m_p.cgs.value

        # reminder: it should be
        #   read from the camera
        bbox_lim = regionsize_kpc / 2.
        bbox = np.array([[-bbox_lim, bbox_lim],
                         [-bbox_lim, bbox_lim],
                         [-bbox_lim, bbox_lim]])

        shape_data = data['density'].shape
        # data should be added with proper units here
        #   (maybe get_units from pymses can be used)

        print data.keys()

        data = dict(density=(mp * data['density'], "g/cm**3"),
                    h2density=(data["density"] * data["H2"], "1/cm**3"),
                    P=(data['P'], "K/cm**3"),
                    P_nt=(data['P_nt'], "K/cm**3"),
                    Z=(data['Z'], ""),
                    velx=(data['velx'], "km/s"),
                    vely=(data['vely'], "km/s"),
                    velz=(data['velz'], "km/s")
                    )

        ds = yt.load_uniform_grid(
            data, shape_data,  length_unit='kpc', bbox=bbox)
    else:
        ds = yt.load_uniform_grid(data, data["density"].shape)

        field = ("h2density")

        def _h2density(field, data):
            try:
                return data["density"] * data["H2"]
            except:
                return data[("stream", "density")] * data[("stream", "H2")]

        # unit is in 1/cc only if convert_unit is properly called when loading
        # in data
        ds.add_field(("stream", "h2density"),
                     function=_h2density, units="g/cm**3")

    dd = ds.all_data()

    if debug:

        prj = yt.ProjectionPlot(ds, 0, 'h2density',
                                center='c', weight_field='density')
        prj.set_unit('h2density', '1/cm**3')
        prj.save('test_h2density_yt_unit_plot.png')
        print 'dump to ', 'test_h2density_yt_unit_plot.png'

        prj = yt.ProjectionPlot(ds, 0, 'density',
                                center='c', weight_field='density')
        #prj.set_unit('density', 'g/cm**3')
        prj.set_unit('density', 'Msun/pc**3')
        #prj.set_unit('density', 'code_mass/code_length**3')
        prj.save('test_density_yt_unit_plot.png')
        print 'dump to ', 'test_density_yt_unit_plot.png'

        print dd['h2density'].max(), dd['h2density'].min()
        print dd['density'].max(), dd['density'].min()

    return ds, dd


def prepare_star_unigrid(data, regionsize_kpc=7., verbose=False, add_unit=False, debug=False):

    import yt
    # prepare unigrid for stars
    bbox_lim = regionsize_kpc / 2.
    bbox = np.array([[-bbox_lim, bbox_lim],
                     [-bbox_lim, bbox_lim],
                     [-bbox_lim, bbox_lim]])

    shape_data = data['velx'].shape

    if add_unit:
        data = dict(mass=(data['mass'], "Msun"),
                    age=(data['epochMyr'], "Myr"),
                    velx=(data['velx'], "km/s"),
                    vely=(data['vely'], "km/s"),
                    velz=(data['velz'], "km/s")
                    )

        ds = yt.load_uniform_grid(
            data, shape_data, length_unit='kpc', bbox=bbox)
    else:

        ds = yt.load_uniform_grid(data, shape_data)

    dd = ds.all_data()

    if verbose:
        print dd['velx'].max(), dd['velx'].min()
        print dd['mass'].max(), dd['mass'].min()
        print dd['epoch'].max(), dd['epoch'].min()
    return ds, dd


def check_hist_h2(data, th_list, ss=None, outdir='./'):
    """

    print summary statistics for H2 density. Plot histogram of H2 density, then plot vertical lines of the cuts.


    Parameters
    ----------
    th_list: list
        list of threshold to make cuts to identify clumps.
    ss: int (optional)
        snapshot number

    """

    if not isinstance(th_list, list):
        th_list = [th_list]

    import matplotlib.pyplot as plt
    aa = (data["density"] * data["H2"]).flatten()
    print np.max(aa), np.min(aa), np.min(aa[aa > 0])
    aa[aa <= 0] = np.min(aa[aa > 0])
    aa = np.log10(aa)

    print np.max(aa), np.min(aa), np.min(aa[aa > 0])

    plt.close('all')
    plt.figure()
    toplot = aa[aa>-7]
    plt.hist(toplot, bins=100) # , normed=True)

    for ele in th_list:
        # plot the ncut positions.
        x = np.log10(ele)
        # plt.plot([x, x], [1, 1.e+7], ls='--', color='k')

    plt.yscale('log')
    # plt.ylabel('Counts')
    plt.xlabel(r"$\log~n_{{\rm H}_2}$ [cm$^{-3}$]")
    plt.minorticks_on()
    plt.tight_layout()
    if ss is None:
        plt.savefig(outdir + 'hist_test.png')
        # plt.savefig(outdir + 'hist_test.pdf')
    else:
        # plt.savefig(outdir + 'hist_test_' + str(ss) + '.png')
        plt.savefig(outdir + 'hist_test_' + str(ss) + '.pdf', bbox_inches='tight')


def check_power(data, size_kpc=7., isnap=28, outdir=''):

    import matplotlib.pyplot as plt

    f_out = outdir + 'test_power.png'

    data_to_process = np.copy(data["density"]) * np.copy(data["H2"])
    n1, n2, n3 = np.shape(data_to_process)

    # comput power
    FT = np.fft.fftn(data_to_process)
    power = FT.real * FT.real + FT.imag * FT.imag

    # set distance
    x = np.linspace(-size_kpc, size_kpc, n1)
    y = np.linspace(-size_kpc, size_kpc, n2)
    z = np.linspace(-size_kpc, size_kpc, n3)
    x, y, z = np.meshgrid(x, y, z)
    dist = np.sqrt(x**2 + y**2 + z**2)

    # 3D -> 1D
    P = power.reshape(np.size(power))
    dist = dist.reshape(np.size(dist))

    intervals = np.linspace(0., size_kpc / 2, n1)

    p = np.histogram(dist, bins=intervals, weights=P)[0]
    pd = np.histogram(dist, bins=intervals)[0]
    pd.astype('float')
    p = p / pd

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(2. * np.pi / intervals[1:], p)
    plt.xscale('log')
    plt.yscale('log')

    ax.set_xlim(0., 2 * np.pi * size_kpc / 2)

    print 'saving to', f_out
    plt.savefig(f_out)

    exit()
