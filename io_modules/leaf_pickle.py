import cPickle as pickle
import sys
sys.path.append('../')

from clump_modules.leafprop import Cloud
from io_modules.manipulate_fetch_gal_fields import get_dx


def load_in_pickled_leaf_singleSS(directory, snapshot_num, pattern="*_10*.p"):
    """

    Load in all the pickled files in a given directory.
    For a single snapshot, load results (leaf clumps) from different ncut.

    Parameters
    ----------
    directory: str
        must be absolute path

    snapshot_num: int
        snapshot number

    Returns
    -------
    ss: dict
        containing all the cloud properties

    """

    import glob
    pfiles = glob.glob(directory + pattern)
    assert len(pfiles) > 0

    ss = {}
    for i, pfnames in enumerate(pfiles):
        starti = pfnames.find(str(snapshot_num) + '/') + 3
        ncut_str = pfnames[starti:starti + pfnames[starti:].find('_')]
        # print ncut_str
        leaf_fields = pickle.load(open(pfnames, "rb"))

        Cloud_dict = {}
        for kkk in leaf_fields.iterkeys():
            Cloud_dict[kkk] = Cloud(get_dx(snapshot_num), leaf_fields[kkk],
                                    int(kkk))
            print Cloud_dict[kkk]
        ss[ncut_str] = Cloud_dict

    return ss


def load_in_pickle_SS(basepath, picklefname, minss, maxss):
    """

    Load in pickled leaf clumps for all snapshots, selected using same criterion (reflected in picklefname).

    Parameters
    ----------
    basepath: str
        must be absolute path

    Returns
    -------
    ss: dict
        containing all the cloud properties, for all snapshots

    """

    ss = {}
    for issnum in range(minss, maxss+1):
        leafdir = basepath + str(issnum) + '/'
        try:
            leaf_fields = pickle.load(open(leafdir + picklefname, "rb"))
        except:
            continue

        Cloud_dict = {}
        for kkk in leaf_fields.iterkeys():
            Cloud_dict[kkk] = Cloud(get_dx(issnum), leaf_fields[kkk], int(kkk))
            print Cloud_dict[kkk]

        ss[str(issnum)] = Cloud_dict

    return ss
