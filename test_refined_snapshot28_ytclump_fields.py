"""

Potential BUG/issue/inconsisntency (may need to update code):  density should be multipled by factor = get_units(ro=ro)['rho'][0]      # 1/cm^3 (not H/cm^3) if we use Fig 6 of Pallottini+17 as n_cut.


http://yt-project.org/doc/analyzing/analysis_modules/clump_finding.html

After getting the finely gridded cube from resample_fields.py, use yt clump finder.


Last mod: 16 July 2018

"""

import numpy as np
import matplotlib.pyplot as plt

import yt
from clump_modules.clump_wrapper import ytclumpfind_H2
from io_modules.manipulate_fetch_gal_fields import import_fetch_gal, prepare_unigrid



if __name__ == '__main__':

    fold_out = 'test_png_28/'
  
    data = import_fetch_gal(isnap = 28)
    ds,dd = prepare_unigrid(data= data)
  
    # -------------- decide on n_cut --------------
    # may be useful to plot the 3D, see test_rendering.py


    # Or see fig 6 Pallottini 2017 for n_H2 cut as starting point, lower right panel, based on Minkowsky function (Euler characteristic).

    # in units of nH2/cc
    n_cut_1 = 10**0.5
    n_cut_2 = 10**-1.5


    # -------------- run clump finder -------------

    master10, leaf10 = ytclumpfind_H2(ds, dd, ("h2density"),
                                      n_cut=0.1,
                                      step=5.0,
                                      N_cell_min=10,
                                      plot=True,
                                      saveplot=True,
                                      fold_out=fold_out)

    # write_clump_index(master10, 0, "master10_clump_hierarchy_H2.txt")
    # write_clumps(master10, 0,  "master10_clumps_H2.txt")

    # import os
    # os.system('cat *_clump_hierarchy_H2.txt')
    # os.system('cat *_clumps_H2.txt')


    # for ind in range(len(leaf10)):
    #     print(leaf10[ind]["h2density"])
    #     print(leaf10[ind].quantities.total_mass())
    #     print(leaf10[ind].quantities.center_of_mass())


    # master20, leaf20 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_1,
    #                                   step=20,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master20, 0, "master20_clump_hierarchy_H2.txt")
    # write_clumps(master20, 0,  "master20_clumps_H2.txt")

    # master30, leaf30 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_1,
    #                                   step=30,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master30, 0, "master30_clump_hierarchy_H2.txt")
    # write_clumps(master30, 0,  "master30_clumps_H2.txt")

    # master70, leaf70 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_1,
    #                                   step=70,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master70, 0, "master70_clump_hierarchy_H2.txt")
    # write_clumps(master70, 0,  "master70_clumps_H2.txt")

    # master100, leaf100 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_1,
    #                                     step=100,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # write_clump_index(master100, 0, "master100_clump_hierarchy_H2.txt")
    # write_clumps(master100, 0,  "master100_clumps_H2.txt")


    # master200, leaf200 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_1,
    #                                     step=200,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # write_clump_index(master200, 0, "master200_clump_hierarchy_H2.txt")
    # write_clumps(master200, 0,  "master200_clumps_H2.txt")


    # # --- repeat for n_cut_2 ---
    # master10, leaf10 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=10,
    #                                   N_cell_min=20,
    #                                   plot=True,
    #                                   saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master10, 0, "master10_clump_hierarchy_H2_cut2.txt")
    # write_clumps(master10, 0,  "master10_clumps_H2_cut2.txt")

    # master20, leaf20 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=20,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master20, 0, "master20_clump_hierarchy_H2_cut2.txt")
    # write_clumps(master20, 0,  "master20_clumps_H2_cut2.txt")

    # master30, leaf30 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=30,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # write_clump_index(master30, 0, "master30_clump_hierarchy_H2_cut2.txt")
    # write_clumps(master30, 0,  "master30_clumps_H2_cut2.txt")

    # master70, leaf70 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                   n_cut=n_cut_2,
    #                                   step=70,
    #                                   N_cell_min=20, plot=True, saveplot=True,
    #                                   fold_out=fold_out)
    # try:
    #     write_clump_index(master70, 0, "master70_clump_hierarchy_H2_cut2.txt")
    #     write_clumps(master70, 0,  "master70_clumps_H2_cut2.txt")
    # except:
    #     pass

    # master100, leaf100 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_2,
    #                                     step=100,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # try:
    #     write_clump_index(master100, 0, "master100_clump_hierarchy_H2_cut2.txt")
    #     write_clumps(master100, 0,  "master100_clumps_H2_cut2.txt")
    # except:
    #     pass

    # master200, leaf200 = ytclumpfind_H2(ds, dd, ("h2density"),
    #                                     n_cut=n_cut_2,
    #                                     step=200,
    #                                     N_cell_min=20, plot=True, saveplot=True,
    #                                     fold_out=fold_out)
    # try:
    #     write_clump_index(master200, 0, "master200_clump_hierarchy_H2_cut2.txt")
    #     write_clumps(master200, 0,  "master200_clumps_H2_cut2.txt")
    # except:
    #     pass

    # # find out leaf_clumps attributes to retreive physical properties
    # aa = leaf10[0]
    # print aa.field
    # print aa.data.fcoords    # grids position
    # print aa.data.icoords    # grids element index
    # print aa.data.icoords[0]
    # print aa.data['h2density'][0]
    # _density = f["rho"].value
    # _h2 = f["H2"].value
    # ii = aa.data.icoords[:, 0]
    # jj = aa.data.icoords[:, 1]
    # kk = aa.data.icoords[:, 2]
    # print (_density[ii, jj, kk] * factor) * _h2[ii, jj, kk]
    # assert round(aa.data['h2density'][0]) == round((_density[ii, jj, kk] * factor * _h2[ii, jj, kk])[0])

