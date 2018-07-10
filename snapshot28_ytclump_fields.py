"""

http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html

After getting the finely gridded cube from resample_fields.py, use yt clump finder.

Last mod: 9 July 2018

NOTE
----
The yt clump finder was initially described in http://adsabs.harvard.edu/abs/2009ApJ...691..441S , but it's changed since then. What it does now is decompose into non-overlapping tiles (stored in a kd-tree), identify contours within a tile, and then connect them across tiles. It does this with an upper and lower bound on a field value, and looking for topologically connected sets.

With yt, connected sets can be identified. This enables the location and analysis of a hierarchy of clumps; in star formation, for instance, this allows the construction of diagrams describing the density at which fragmentation occurs.

- extract_connected_sets() differ slightly from Clump()
    - extract_connected_sets(): get topologically connected sets. These sets are identified by examining cells between two threshold values and connecting them.
    - Clump(): uses a contouring algorithm to identified topologically disconnected structures within a dataset. This works by first creating a single contour over the full range of the contouring field, then continually increasing the lower value of the contour until it reaches the maximum value of the field. As disconnected structures are identified as separate contours, the routine continues recursively through each object, creating a hierarchy of clumps. Individual clumps can be kept or removed from the hierarchy based on the result of user-specified functions, such as checking for gravitational boundedness.

"""

import numpy as np
import matplotlib.pyplot as plt

import yt
from yt.analysis_modules.level_sets.api import Clump, find_clumps, get_lowest_clumps, write_clump_index, write_clumps
import h5py

Plot_stuff = False
debug = False

fold_out   = 'test_png'


# convert from code unit density to g/cc (depending on how fetch_gal.py is
# implemented.)
convert_unit = True

if debug:
    namefile = "output/output_00028/info_00028.txt"
    myfile = namefile

    print "Loading file,", myfile
    pf = load(myfile)


f = h5py.File("snapshot28_center_fields0123456-15_resampled.h5", "r")
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


if convert_unit:
    import pymses
    from pymses.utils import constants as C

    ro = pymses.RamsesOutput("output", 28)

    from fetch_gal_fields import get_units
    factor_density = get_units(ro=ro)['rho'][0]      # 1/cm^3 (not H/cm^3)
    density *= factor_density
    print density.max()

    factor_vel = get_units(ro=ro)['vel'][0]
    velx *= factor_vel
    vely *= factor_vel
    velz *= factor_vel
    print velx.max(), vely.max(), velz.max()

    factor_P = get_units(ro=ro)['P'][0]
    Pressure *= factor_P
    P_nt *= factor_P
    print np.log10(Pressure.max()), np.log10(P_nt.max())


data = dict(density=density, H2=H2,
            P=Pressure,
            P_nt=P_nt,
            Z=metallicity,
            velx=velx,
            vely=vely,
            velz=velz)

ds = yt.load_uniform_grid(data, f["rho"].shape)
dd = ds.all_data()

# make a derived field, call h2density (for yt Clump() to work)
def _h2density(field, data):
    try:
        return data["density"] * data["H2"]
    except:
        return data[("stream", "density")] * data[("stream", "H2")]


from yt.units import dimensions
ds.add_field(("stream", "h2density"), function=_h2density, units="code_mass/code_length**3")
print dd['h2density'].max()

assert (dd['H2'] * dd['density']).max() == dd['h2density'].max()


def ytclumpfind_H2(ds, dd, field, n_cut, step=10, N_cell_min=20, save=False, plot=True, saveplot=None, fold_out='./'):
    '''

    The way it's implemented now only works for single "density" field.

    Parameter
    ---------
    ds: yt StreamDataset

    dd: YTRegion

    field: tuple
        tuple of str, e.g., ("gas", "density") or ("io", "....") or ("gas", "averaged_density"), etc... or just ("density") or ("stream", "density")

    n_cut: float or int
        defines lowest contour level to start searching from.
        e.g., density to cut the clouds.

    step: int
        multiplicative interval between subsequent contours

    N_cell_min: int
        min. number of cell s.t. clump identified is not spurious

    save: Boolean
        if true, save the clump tree as a reloadable dataset

    plot: Boolean
        if true, will plot leaf clump

    saveplot: boolean
        if true, will save figure instead of showing it

    fold_out: str
        directory to save figures


    Return
    ------
    master_clump: the top of a hierarchy of clumps
    leaf_clumps: list of individual clumps that have no children of their own


    '''

    if plot:
        assert saveplot is not None

    # c_min = 10**np.floor(np.log10(dd[field]).min()  )
    # c_max = 10**np.floor(np.log10(dd[field]).max()+1)
    if n_cut < 1.e-5:
        n_cut = 1.0    # to make sure whatever comes out after multiplicative by step won't be too small
    c_min = n_cut
    c_max = (dd[field]).max()

    # assert np.isnan(c_min) is False and np.isnan(c_max) is False

    print "min/max value for finding contours: ", c_min, c_max

    # this "base clump" just  covers the whole domain.
    master_clump = Clump(dd, field)
    # weed out clumps < N_cell_min cells.
    master_clump.add_validator("min_cells", N_cell_min)

    find_clumps(master_clump, c_min, c_max, step)

    if save:

        fn = master_clump.save_as_dataset(
            fields=list(field)),  # "particle_mass"])
        # # To reload the clump dataset
        # cds = yt.load(fn)
        # leaf_clumps_reloaded = cds.leaves

    # traverse clump hierarchy to get list of all 'leaf' clumps, which are the
    # individual clumps that have no children of their own
    leaf_clumps = get_lowest_clumps(master_clump)

    def plotclumps(ds, field=field, saveplot=saveplot, fold_out=fold_out):
        """ overplot the clumps found (specifically the leaf_clumps) along 3 images, each created by projecting onto x-, y-, and z-axis. """

        axes = {'0': 'x', '1': 'y', '2': 'z'}

        for kk, vv in axes.iteritems():

            prj = yt.ProjectionPlot(ds,
                                    int(kk),
                                    field,
                                    center='c')
            prj.annotate_clumps(leaf_clumps)
            if saveplot:
                prj.save(fold_out + 'clumps1_' + str(int(n_cut)) + '_' + str(int(step)) + '-' + str(int(N_cell_min)) + '_' + vv + 'axis.png')
            else:
                prj.show()

    if plot:
        plotclumps(ds, saveplot=saveplot, fold_out=fold_out)

    return master_clump, leaf_clumps


if __name__ == '__main__':
    # -------------- decide on n_cut --------------
    # fig 6 Pallottini 2017 for n_H2 cut as starting point, lower right panel, based on Minkowsky function (Euler characteristic).

    # in units of nH2/cc
    n_cut_1 = 10**0.5
    n_cut_2 = 10**-1.5

    # -------------- run clump finder -------------
    master5, leaf5 = ytclumpfind_H2(ds, dd, ("h2density"),
                                    n_cut=n_cut_1,
                                    step=5,
                                    N_cell_min=3,
                                    plot=True,
                                    saveplot=True,
                                    fold_out=fold_out)

    master10, leaf10 = ytclumpfind_H2(ds, dd, ("h2density"),
                                      n_cut=n_cut_1,
                                      step=10,
                                      N_cell_min=3,
                                      plot=True,
                                      saveplot=True,
                                      fold_out=fold_out)

    print(master10.children)
    print(master10.children[0]['h2density'])    # children
    # print(master10.children[0].children[0]['h2density'])   # grand-children; note not necessary that children of master has a grandchild


    for ind in range(len(leaf10)):
        print(leaf10[ind]["h2density"])
        print(leaf10[ind].quantities.total_mass())
        print(leaf10[ind].quantities.center_of_mass())


    # --- repeat for n_cut_2 ---
    master5, leaf5 = ytclumpfind_H2(ds, dd, ("h2density"),
                                    n_cut=n_cut_2,
                                    step=5,
                                    N_cell_min=3,
                                    plot=True,
                                    saveplot=True,
                                    fold_out=fold_out)

    master10, leaf10 = ytclumpfind_H2(ds, dd, ("h2density"),
                                      n_cut=n_cut_2,
                                      step=10,
                                      N_cell_min=3,
                                      plot=True,
                                      saveplot=True,
                                      fold_out=fold_out)

    # to retreive physical properties, i.e., need to get cell indice of clumps (still need some more work here..)
    aa = leaf10[0]
    print aa.field
    ii = aa.data.icoords[:, 0]
    jj = aa.data.icoords[:, 1]
    kk = aa.data.icoords[:, 2]

    print aa.data.icoords[0]
    print aa.data['h2density'][0]

    density_aa = density[ii, jj, kk]
    H2_aa = H2[ii, jj, kk]
    Pressure_aa = Pressure[ii, jj, kk]
    P_nt_aa = P_nt[ii, jj, kk]
    metallicity = metallicity[ii, jj, kk]
    velx = velx[ii, jj, kk]
    vely = vely[ii, jj, kk]
    velz = velz[ii, jj, kk]

    _h2_aa = aa.data['h2density']
    h2_aa = (density_aa * factor_density) * H2_aa
    assert round(_h2_aa[0]) == round(h2_aa[0])

    import pdb; pdb.set_trace()




# then... look at fields in these corresponding grid points for the "clumps" identified.....


# seems to need lower steps to identify more structures.. will do after cleaning up code maybe, or write in another python script for this.