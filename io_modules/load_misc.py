import numpy as np

def load_SFR_perSS(filename="/mnt/home/daisyleung/mc_eor/star_particles/SFRs.txt"):
    """ load in globally-integrated SFR calculated for each snapshot, from a file saved before.
    """
    sfr = {}

    ss, sfr_4myr = np.loadtxt(filename, unpack=True, usecols=(0, 2))

    for idx, iii in enumerate(ss):
        sfr[str(int(iii))] = sfr_4myr[idx]

    return sfr

