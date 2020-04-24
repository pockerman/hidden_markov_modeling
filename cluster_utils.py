import numpy as np
from sklearn.neighbors import KernelDensity
from helpers import Error


def build_cluster_densities(clusters, windows, **kwargs):
  """
  Establish the probability destributions underlying
  the data for each cluster

  Parameters
  ----------
  clusters : TYPE
    DESCRIPTION.
  **kwargs : TYPE
    DESCRIPTION.

  Returns
  -------
  None.

  """
  if kwargs["cluster_distribution"]["name"] == "kde":

    for cluster in clusters:

      indeces = cluster.indexes

      arr = np.empty((0, 2), float)

      for idx in indeces:
        window = windows[idx]
        mu1, mu2 = window.get_rd_stats(statistics="mean", name="both")
        arr = np.append(arr, np.array([[mu1, mu2]]), axis=0)

      kde = KernelDensity(kernel=kwargs["config"]["kernel"],
                          bandwidth=kwargs["config"]["bandwidth"])
      kde.fit(arr)
      cluster.density = kde



  raise Error("Invalid cluster distribution method")


def clusters_statistics(clusters, windows):
  """
  Claculate various statistics for the windows
  clustered in clusters

  Parameters
  ----------
  clusters : list of lists
    Contains the clustered window indexes. There are
    len(clusters) clusters
  windows : list of Window objects
    DESCRIPTION.

  Returns
  -------
  a map with the calculated statistics for each
  cluster

  """

  statistics = {}

  for c in range(len(clusters)):

    statistics[c] = clusters[c].get_statistics(windows=windows, statistic="all")
  return statistics

