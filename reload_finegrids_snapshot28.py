"""

reload into yt object after we have regridding the AMR data to uniform grids.

Last mod: 5 July 2018


"""

import yt
import h5py

f = h5py.File("snapshot28_center_densityfield_resampled.h5", "r")
density = f["density"].value


# convert from code unit density to g/cc
convert_unit = True

if convert_unit:
    import pymses
    from pymses.utils import constants as C

    ro = pymses.RamsesOutput("output", 28)
    factor = ro.info["unit_density"].express(C.H_cc)
    density *= factor
    print density.max()

data = dict(density=density)
ds = yt.load_uniform_grid(data, f["density"].shape)

dd = ds.all_data()
# dd['density']