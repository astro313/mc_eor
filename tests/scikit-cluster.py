'''

Scikit-learn clustering methods.

'''

print(__doc__)
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
import h5py

f = h5py.File("snapshot28_center_densityfield_resampled.h5", "r")
density = f["density"].value

# convert from code unit density to g/cc (depending on how fetch_gal.py is implemented.)
convert_unit = True

if convert_unit:
    import pymses
    from pymses.utils import constants as C

    ro = pymses.RamsesOutput("output", 28)
    factor = ro.info["unit_density"].express(C.H_cc)
    density *= factor
    print density.max()


'''

DBSCAN: views clusters as areas of high density separated by areas of low density.
-  uses ball trees and kd-trees to determine the neighborhood of points
- clusters found by DBSCAN can be any shape, as opposed to k-means which assumes that clusters are convex shaped.
- Higher min_samples or lower eps indicate higher density necessary to form a cluster.
    -  because a core sample = a sample in the dataset s.t. there exist "min_samples" number of other samples within a distance of "eps".
- A cluster is a set of core samples, can be built recursively by taking a core sample, finding all of its neighbors that are core samples, finding all of their neighbors that are core samples, and so on.

- Any core sample is part of a cluster, by definition. Any sample that is not a core sample, and is > eps away from any core sample = outlier.


'''

db = DBSCAN(eps=0.3, min_samples=10).fit(density)
labels = db.labels_
core_samples_mask = np.zeros_like(labels, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of clusters: %d' % n_clusters_)
print("Silhouette Coefficient: %0.3f"
      % metrics.silhouette_score(density, labels))


# -------- Plot result -----
# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()

# ----------