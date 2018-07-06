'''

resample AMR points to finest resolution


last mod: 5 July 2018

'''

import numpy as np

ds = np.load('snapshot28_center.npz')    # saved in fetch_gal.py
dx_vector = ds['dx_vector']
density_vector = ds['density_vector']
loc_vector = ds['loc_vector']

originalLevels = np.log2(1./np.unique(dx_vector))
numLevel = len(originalLevels)   # or numLevel = len(np.unique(dx_vector))
levels = originalLevels               # let's just keep as
originalLevels_vector = np.log2(1./dx_vector)
(originalLevels_vector == np.sort(originalLevels_vector)).all() == True


# ii = np.argsort(-dx_vector)
# assert dx_vector[ii[0]] > dx_vector[ii[-1]]    # make sure ii starts w/ coarsest resolution


levels_ii = originalLevels_vector
levels_ii = np.array(levels_ii, dtype=int)
levels = np.array(levels, dtype=int)

d = {}
for LL in levels:
    # find indices where they correspond to each of the levels
    mask = levels_ii == LL
    d[str(LL)] = mask
    print LL, np.max(levels_ii[mask]), np.min(levels_ii[mask])


# check originalLevels_vector[ii] is properly saved to the right dict key
ilowestLevel = d[str(int(levels.min()))]
# assert (np.sort(ilowestLevel) == np.sort(np.where(originalLevels_vector == levels.min())[0])).all() == True

# coarset level
Llow = str(int(levels.min()))
pos  = loc_vector[d[Llow],:]
II   = density_vector[d[Llow]]

ilev = levels_ii[d[Llow]]
assert ilev.min() == ilev.max() == int(Llow)

print density_vector
print II


originalSize = 0.0015
highestRes = 2.**(levels.max() * -1)
N = originalSize/highestRes
print N

# make sure it's even number
# N = 196
N = int(N)

Ninit = N/len(levels)/2
density_cube = np.zeros((Ninit, Ninit, Ninit))

imatrix = np.ones((2, 2, 2))

for lll in np.sort(levels):

    print lll
    density_cube = np.kron(density_cube, imatrix)

    pos = loc_vector[d[str(lll)], :]

    xx = (pos[:, 0] - loc_vector[:, 0].min())/(loc_vector[:, 0].max() - loc_vector[:, 0].min())
    yy = (pos[:, 1] - loc_vector[:, 1].min())/(loc_vector[:, 1].max() - loc_vector[:, 1].min())
    zz = (pos[:, 2] - loc_vector[:, 2].min())/(loc_vector[:, 2].max() - loc_vector[:, 2].min())

    print '  ',xx.min(), xx.max()
    print '  ',yy.min(), yy.max()
    print '  ',zz.min(), zz.max()
    #import pdb; pdb.set_trace()

    xx   = xx * density_cube.shape[0]
    yy   = yy * density_cube.shape[1]
    zz   = zz * density_cube.shape[2]

    print '  ',xx.min(), xx.max(), yy.min(), yy.max(), zz.min(), zz.max()

    xpos = np.array(xx, dtype=int)
    ypos = np.array(yy, dtype=int)
    zpos = np.array(zz, dtype=int)

    xpos[xpos == density_cube.shape[0]] = density_cube.shape[0] -1
    ypos[ypos == density_cube.shape[1]] = density_cube.shape[1] -1
    zpos[zpos == density_cube.shape[2]] = density_cube.shape[2] -1
    density_cube[xpos, ypos, zpos] = density_vector[d[str(lll)]]


size = density_cube.shape[0]

import matplotlib.pyplot as plt
plt.imshow(np.log10(density_cube[:, :, 2**4:].sum(axis=0)))   # 2**4: to chop off spurious edges
plt.show()


np.savez_compressed('snapshot28_center_densityfield_resampled',
                    density_cube=density_cube,
                    highestRes=highestRes)
