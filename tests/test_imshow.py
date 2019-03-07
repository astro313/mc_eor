from yt.visualization.fixed_resolution import FixedResolutionBuffer
import h5py
import yt
import numpy as np

snapshot_num = 28
f = h5py.File("snapshot" + str(snapshot_num) + "_center_fields0123456-15_resampled.h5", "r")
# careful sometimes i used "density" (e.g., resample.py), see
# resample_fields.py to make sure
density     = f["rho"].value
H2          = f["H2"].value
Pressure    = f["P"].value
P_nt        = f["P_nt"].value
metallicity = f["Z"].value
velx        = f["vel_x"].value
vely        = f["vel_y"].value
velz        = f["vel_z"].value


import pymses
from pymses.utils import constants as C

ro = pymses.RamsesOutput("output", snapshot_num)

from io_modules.manipulate_fetch_gal_fields import get_units
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


data = dict(density = density, H2=H2,
            P       = Pressure,
            P_nt    = P_nt,
            Z       = metallicity,
            velx    = velx,
            vely    = vely,
            velz    = velz
            )

ds = yt.load_uniform_grid(data, f["rho"].shape)
dd = ds.all_data()

# make a derived field, call h2density (for yt Clump() to work)
def _h2density(field, data):
    try:
        return data["density"] * data["H2"]
    except:
        return data[("stream", "density")] * data[("stream", "H2")]


from yt.units import dimensions
ds.add_field(("stream", "h2density"), function=_h2density,
             units="code_mass/code_length**3")
print dd['h2density'].max()

assert (dd['H2'] * dd['density']).max() == dd['h2density'].max()

# prj = yt.ProjectionPlot(ds,
#                         0,
#                         "h2density",
#                         center='c')
prj = ds.proj("h2density", 0)
# prj.annotate_clumps(leaf_clumps)

frb = FixedResolutionBuffer(prj, (-1,1,-1,1), (800,800))

my_image = frb["h2density"]
fits.writeto("my_images.fits", my_image)



