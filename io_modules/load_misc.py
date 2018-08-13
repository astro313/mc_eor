import numpy as np
import cPickle as pickle

import os
here     = os.path.dirname(os.path.abspath(__file__))
datapath = here[:here.rfind('/')]+'/star_particles/'

def load_SFR_perSS(filename=datapath+"SFRs.txt"):
    """ load in globally-integrated SFR calculated for each snapshot, from a file saved before.
    """
    sfr = {}

    ss, sfr_4myr = np.loadtxt(filename, unpack=True, usecols=(0, 2))

    try:
      for idx, iii in enumerate(ss):
        sfr[str(int(iii))] = sfr_4myr[idx]
    except TypeError:
      sfr[str(int(ss))] = sfr_4myr

    return sfr



def get_camera_from_file(f_camera):

    """ 

    Parameters
    ----------
    f_camera: str
       

    Returns
    -------
    camera: dict
       containing at least 'center' and 'region_size' for snapshots 

   """ 

    with open(f_camera, 'rb') as f:
        camera = pickle.load(f)

    return camera
