import numpy as np
from sklearn.neighbors import KernelDensity
from helpers import Error
from helpers import INFO


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
  print("{0} Type of cluster density fitted: {1}".format(INFO,
                                                         kwargs["name"]))
  if kwargs["name"] == "kde":

    print("{0} Type of cluster density kernel: {1}".format(INFO,
                                                           kwargs["config"]["kernel"]))

    for cluster in clusters:

      arr = _form_cluster_2d_array(cluster=cluster, windows=windows)
      kde = KernelDensity(kernel=kwargs["config"]["kernel"],
                          bandwidth=kwargs["config"]["bandwidth"])
      kde.fit(arr)
      cluster.density = kde
    return clusters

  raise Error("Invalid cluster distribution method")


def save_clusters_desnity(clusters, windows, **kwargs):

  for cluster in clusters:
    filename = "cluster_" + str(cluster.cidx) + "_density.txt"
    save_cluster_density(cluster=cluster, windows=windows,
                         filename=filename, **kwargs)


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


def save_cluster_density(cluster, windows, filename, kwargs):


    if kwargs["name"] == "kde":

      arr = _form_cluster_2d_array(cluster=cluster, windows=windows)
      log_probs = cluster.density.score_logs(arr)

      with open(filename, 'w') as file:
        file.write(str(log_probs))
      return

    raise Error("Invalid cluster distribution method")


def _form_cluster_2d_array(cluster, windows):
  indeces = cluster.indexes
  arr = np.empty((0, 2), float)

  for idx in indeces:
    window = windows[idx]
    mu1, mu2 = window.get_rd_stats(statistics="mean", name="both")
    arr = np.append(arr, np.array([[mu1, mu2]]), axis=0)
  return arr


