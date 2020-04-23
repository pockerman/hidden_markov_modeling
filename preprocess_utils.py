"""
Preprocessing utilities
"""

from pomegranate import*
from scipy import stats
import numpy as np

from exceptions import Error

from helpers import listify_dicts_property
from helpers import WindowState
from helpers import flat_windows
from helpers import INFO

VALID_DISTS = ['normal', 'uniform',
               'poisson','discrete',]


def fit_distribution(data, dist_name="normal", **kwargs):

    """
    Fits a distribution within the given dataset
    """
    if dist_name not in VALID_DISTS:
        raise Error("Invalid distribution name. \
                    Name '{0}' not in {1}".format(dist_name, VALID_DISTS))

    if dist_name == 'normal':
        dist = NormalDistribution.from_samples(data)
        return dist
    elif dist_name == 'uniform':
        dist = UniformDistribution.from_samples(data)
        return dist
    elif dist_name == 'poisson':
        dist = PoissonDistribution.from_samples(data)
        return dist
    elif dist_name == 'discrete':
        dist = DiscreteDistribution.from_samples(data)
        return dist

def compute_statistic(data, statistics):

  if statistics == "mean":
    return np.mean(data)
  elif statistics == "var":
    return np.var(data)
  elif statistics == "median":
    return np.median(data)
  elif statistics == "min":
    return np.amin(data)
  elif statistics == "max":
    return np.amax(data)
  elif statistics == "mode":
    return stats.mode(data, axis=None).mode[0]
  elif statistics == "q75":
    return np.percentile(data, [75])
  elif statistics == "q25":
    return np.percentile(data, [25])
  elif statistics == "q50":
    return np.percentile(data, [50])
  elif statistics == "all":

          mean = np.mean(data)
          var = np.var(data)
          median = np.median(data)
          min_ = np.amin(data)
          max_ = np.amax(data)
          mode = stats.mode(data, axis=None).mode[0]
          q75, q50, q25 = np.percentile(data, [75, 50, 25])
          return {"mean": mean, "var": var,
                  "median": median,
                  "min": min_,
                  "max": max_,
                  "mode": mode,
                  "iqr": q75 - q25,
                  "q75": q75,
                  "q25": q25,
                  "q50": q50}

def zscore_outlier_removal(windows, config):

  newwindows = []

  for window in windows:
    mu = window.get_rd_stats(statistics="mean", name="both")

    sigma_wga = np.sqrt(config["statistics"]["wga_w"]["var"])
    zscore_wga = (mu[0] - config["statistics"]["wga_w"]["mean"])/sigma_wga

    sigma_no_wga = np.sqrt(config["statistics"]["n_wga_w"]["var"])
    zscore_no_wga = (mu[1] - config["statistics"]["n_wga_w"]["mean"])/sigma_no_wga


    if zscore_wga < - config["sigma_factor"] or\
      zscore_wga > config["sigma_factor"]:
        continue
    elif zscore_no_wga < - config["sigma_factor"] or\
      zscore_no_wga > config["sigma_factor"]:
        continue
    else:
      newwindows.append(window)

  return newwindows


def remove_outliers(windows, removemethod, config):

  if removemethod == "zscore":
    return zscore_outlier_removal(windows=windows, config=config)

  raise Error("Unknown outlier removal method: {0}".format(removemethod))


def build_clusterer(data, nclusters, method, **kwargs):

  """
  A simple wrapper to various clustering approaches.
  Cluster the given data into nclusters by using the
  specified method. Depending on the specified method
  different packages may be required and different
  argumens are expected in the kwargs dict.
  """

  features = kwargs["clusterer"]["config"]["features"]
  windows = []

  print("{0} cluster features used {1}".format(INFO, features))

  for window in data:

    window_data = window.get_rd_stats(statistics="all")
    window_values =[]

    for feature in features:
      window_values.append(window_data[0][feature])
      window_values.append(window_data[1][feature])

    windows.append(window_values)

  if method == "kmeans":

    from sklearn.cluster import KMeans
    clusterer = KMeans(n_clusters=nclusters)

    clusterer.fit(windows)
    return clusterer
  elif method == "kmedoids":
    from pyclustering.cluster.kmedoids import kmedoids
    from pyclustering.utils.metric import type_metric
    from pyclustering.utils.metric import  distance_metric

    if kwargs["clusterer"]["config"]["metric"] == "MANHATAN":
      t_metric= type_metric.MANHATTAN
    elif kwargs["clusterer"]["config"]["metric"] == "EUCLIDEAN":
      t_metric = type_metric.EUCLIDEAN
    else:
      raise Error("Invalid metric specified %s"%kwargs["clusterer"]["config"]["metric"])

    metric = distance_metric(metric_type=t_metric)


    initial_index_medoids=[]
    if kwargs["clusterer"]["config"]["init_cluster_idx"] == "random_from_data":
      import random

      for c in range(nclusters):
        idx = random.randint(0, len(windows)-1)

        if idx in initial_index_medoids:

          # try ten times before quiting
          for time in range(10):
            idx = random.randint(0, len(windows)-1)

            if idx in initial_index_medoids:
              continue
            else:
              initial_index_medoids.append(idx)
              break

        else:
          initial_index_medoids.append(idx)
    else:
        initial_index_medoids=kwargs["clusterer"]["config"]["init_cluster_idx"]


    clusterer  = kmedoids(data=windows,
                          initial_index_medoids=initial_index_medoids,
                          metric=metric)
    clusterer.process()
    return clusterer, initial_index_medoids


  raise Error("Invalid clustering method: " + method )

