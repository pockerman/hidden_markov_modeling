import numpy as np
from sklearn.neighbors import KernelDensity
from pomegranate import *
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

      arr = _form_cluster_2d_array(cluster=cluster)
      kde = KernelDensity(kernel=kwargs["config"]["kernel"],
                          bandwidth=kwargs["config"]["bandwidth"])
      kde.fit(arr)
      cluster.density = kde
  elif kwargs["name"] == "gmm":

    # distributions for the first component
    wga_dist = kwargs["config"]["wga_dist"]
    wga_weights = kwargs["config"]["wga_weights"]

    print("{0} Distributions for WGA {1} ".format(INFO, wga_dist))

    non_wga_dist = kwargs["config"]["no_wga_dist"]
    no_wga_weights = kwargs["config"]["no_wga_dist_weights"]

    print("{0} Distributions for NON WGA {1} ".format(INFO, wga_dist))

    for cluster in clusters:
      indeces = cluster.indexes

      wga_data = np.empty((1,0), float)
      no_wga_data = np.empty((1,0), float)
      windows = clusterwindows

      for idx in indeces:
        window = windows[idx]
        mu1, mu2 = window.get_rd_stats(statistics="mean", name="both")
        wga_data = np.append(wga_data, np.array(mu1))
        no_wga_data = np.append(no_wga_data,
                               np.array(mu2))


      # collected the data create the GMM for each
      # component in the cluster
      params={"mean": np.mean(wga_data),
             "std": np.std(wga_data),
             "std_factor" : 3}
      wga_gmm = \
        GeneralMixtureModel(_get_distributions_list_from_names(wga_dist, params),
                            weights=wga_weights)

      cluster.wga_density = wga_gmm

      params={"mean": np.mean(no_wga_data),
             "std": np.std(no_wga_data),
             "std_factor" : 3}

      non_wga_density = \
        GeneralMixtureModel(_get_distributions_list_from_names(non_wga_dist),
                            weights=no_wga_weights )

      cluster.no_wga_density = non_wga_density

  else:
    raise Error("Invalid cluster distribution method")

  return clusters


def save_clusters_desnity(clusters, **kwargs):

  for cluster in clusters:
    filename = "cluster_" + str(cluster.cidx) + "_density.txt"
    save_cluster_density(cluster=cluster, filename=filename, **kwargs)


def clusters_statistics(clusters):
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

    statistics[c] = clusters[c].get_statistics(statistic="all")
  return statistics


def save_cluster_density(cluster, filename, **kwargs):


    if kwargs["name"] == "kde":

      arr = _form_cluster_2d_array(cluster=cluster)
      log_probs = cluster.density.score_samples(arr)

      with open(filename, 'w') as file:
        file.write(str(log_probs))
      return

    raise Error("Invalid cluster distribution method")


def _form_cluster_2d_array(cluster):
  indeces = cluster.indexes
  arr = np.empty((0, 2), float)
  windows = cluster.windows
  for idx in indeces:
    window = windows[idx]
    mu1, mu2 = window.get_rd_stats(statistics="mean", name="both")
    arr = np.append(arr, np.array([[mu1, mu2]]), axis=0)
  return arr


def _get_distributions_list_from_names(dists_name, params):

  dists = []

  for name in dists_name:
    if name == "normal":
      dists.append(NormalDistribution(params["mean"], params["std"]))
    elif name == "poisson":
      dists.append(PoissonDistribution(params["mean"]))
    elif name == "uniform":
      dists.append(UniformDistribution(params["mean"] - params["std_factor"]*params["std"],
                                       params["mean"] + params["std_factor"]*params["std"]))
    else:
      raise Error("Name {0} is an unknown distribution ".format(name))
  return dists

