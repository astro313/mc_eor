'''

last mod: 16 July 2018

'''

import numpy as np
import pymses

def get_dx(ssnum):
    # hard-code to get dx for now...
    # saved in fetch_gal_fields.py
    ro = pymses.RamsesOutput("output", ssnum)
    ds = np.load('snapshot'+ str(ssnum) +'_center_fields0123456-15.npz')
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

def import_fetch_gal(isnap=28, folder_ramsesdata='output', tag_h5file="_center_fields0123456-15_resampled.h5", verbose=True, convert = True):

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


def import_fetch_stars(isnap=28, folder_ramsesdata='output', tag_h5file="_center_stars_resampled.h5", \
                       verbose=True, convert=True):

    """

    Returns
    -------
    dataDict: dict
        containing star fields (mass and epoch) in converted units

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
    mass = f["mass"].value
    epoch = f["epoch"].value

    from io_modules.pymses_helper import calculate_age_stars
    UniverseAgeStarFormationMyr = calculate_age_stars(ro_in=ro, dset_in ={'epoch': epoch})   # input for calculate_age_stars() should be in code unit

    mask = mass>0
    epoch[mask] = UniverseAgeStarFormationMyr[mask]
    
    if verbose:
        print 'max and min of Universe age when stars are formed [Gyr] '
        print epoch.max()/1.e3, epoch[epoch>0].min()/1.e3 

    if convert:
      factor_mass, unit_mass = get_units(ro=ro)['mass']      
      mass *= factor_mass
      if(verbose):
          print 'max mass in log10'
          print np.log10(mass.max()), unit_mass
          print ''
          print 'sum of all masses in log: Msun'
          print np.log10(mass.sum())

    dataDict = dict(mass=mass, 
                epoch=epoch)
    if True:
        if verbose:
            print 'Clipping variables'
        dataDict["mass"][dataDict["mass"] <= 0] = np.min(dataDict["mass"][dataDict["mass"] > 0])

    if verbose:
        for var in ['mass']:
            print ' Mass range in log10, after clipping: ', var, np.log10(np.max(dataDict[var])), np.log10(np.min(dataDict[var]))

    return dataDict


def prepare_unigrid(data, verbose=False, add_unit= False):

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

    if( add_unit):
      print 'Feature not yet fully implemented'
      # see http://yt-project.org/doc/examining/generic_array_data.html
      # for reference
      
      unit_base = {
          'UnitLength_in_cm'         : 3.08568e+21
         ,'UnitMass_in_g'            : 1.6726219e-24,
          }#,
                 # 'UnitVelocity_in_cm_per_s' :      100000}
      # reminder: it should be
      #   read from the camera
      bbox_lim = 7. # kpc
      bbox_lim = bbox_lim/2.
      bbox     = np.array([[-bbox_lim, bbox_lim],
                           [-bbox_lim, bbox_lim],
                           [-bbox_lim, bbox_lim]])

      data = dict(density = (data['density'], "g/cm**3"))

      ds = yt.load_uniform_grid(data, data["density"].shape,  length_unit='kpc', bbox=bbox)
    else:
      ds = yt.load_uniform_grid(data, data["density"].shape )

    dd = ds.all_data()
    ds.add_field(("stream", "h2density"), function=_h2density, units="g/cm**3")  # unit is in 1/cc only if convert_unit is properly called when loading in data
    assert (dd['H2'] * dd['density']).max() == dd['h2density'].max()

    if add_units:
        prj = yt.ProjectionPlot(ds, 0, field_select,
                    center='c', weight_field='h2density')
        prj.save('test_h2density_yt_unit_plot.png')
        import sys; sys.exit('Exiting..')

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
    plt.hist(aa, bins=100)
    for ele in th_list:
        x = np.log10(ele)
        plt.plot([x, x], [1, 1.e+7], ls='--', color='k')
    plt.yscale('log')
    if ss is None:
        plt.savefig(outdir + 'hist_test.png')
    else:
        plt.savefig(outdir + 'hist_test_' + str(ss) + '.png')
