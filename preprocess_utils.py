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

VALID_DISTS = ['normal', 'uniform', 'poisson',
               'discrete',]


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
    mu = window.get_rd_stats(statistics="mean")

    sigma = np.sqrt(config["statistics"]["var"])
    zscore = (mu - config["statistics"]["mean"])/sigma



    if zscore < - config["sigma_factor"] or\
      zscore > config["sigma_factor"]:
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

  if method == "kmeans":

    from sklearn.cluster import KMeans
    clusterer = KMeans(n_clusters=nclusters)
    windows = flat_windows(data)
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

    windows = flat_windows(data)

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


    clusterer = kmedoids(data=windows,
                         initial_index_medoids=initial_index_medoids,
                         metric=metric)
    clusterer.process()
    return clusterer
  elif method == "wmode":
    return mode_window_clusterer(windows=data,
                                 normal_rd=kwargs["normal_rd"],
                                 delete_rd=kwargs["delete_rd"],
                                 insert_rd=kwargs["insert_rd"])


  raise Error("Invalid clustering method: " + method )


def mode_window_clusterer(windows, normal_rd,
                          delete_rd, insert_rd):

  for window in windows:
    mode = window.get_rd_stats(statistics="mode")

    if mode == normal_rd:
      window.set_state(WindowState.NORMAL)
    elif mode == delete_rd:
      window.set_state(WindowState.DELETE)
    elif mode == insert_rd:
      window.set_state(WindowState.INSERT)
  return windows


"""
def calculate_windows_zscore(windows):
  '''
  Returns an array of Z-scores for every window
  in the given windows list

  Parameters
  ----------
  windows : list of Window instances
    DESCRIPTION.

  Returns
  -------
  list of floats representing the z-score for each window

  '''

  mu = windows_rd_statistics(windows, statistic="mean")
  var = windows_rd_statistics(windows, statistic="var")

  z_scores = []

  for window in windows:
    zscore = (window.get_rd_stats(statistics="mean") - mu)/np.sqrt(var)
    z_scores.append(zscore)

  return z_scores

def windows_tails_p(zscores, interval_length):
  '''
  Calculate windows upper and lower tail probabilities arranged
  into window intervals of length interval_length

  Parameters
  ----------
  zscores :  list
    DESCRIPTION. The score corresponding to each window
  interval_length : int
    Length of window interval

  Returns
  -------
  upper_ps, lower_ps : list of upper probabilities,
                        list of lower probabilities

  '''

  if interval_length > len(zscores):
    raise Error(("You request an interval length {0} not supported by"+
                " the total number of windows {1}").format(interval_length, len(zscores)))

  local_probs_up = []
  local_probs_low = []
  upper_ps = []
  lower_ps = []

  def get_upper_lower_prob_from_score(score):

      # calculate the probability density
      # up to the zscore point
      # i.e caluclate P(Z <= z_i)
      prob = stats.norm.cdf(score)

      # p(Z > z_i) = 1.0 - p(Z <= z_i)
      return 1.0 - prob, prob

  if interval_length == 1:
    for i, zscore in enumerate(zscores):

      # calculate the probability density
      # up to the zscore point
      # i.e caluclate P(Z <= z_i)
      prob = stats.norm.cdf(zscore)

      up, low = get_upper_lower_prob_from_score(score=zscore)

      # p(Z > z_i) = 1.0 - p(Z <= z_i)
      local_probs_up.append({"idx":i,
                             "prob": up})

      local_probs_low.append({"idx":i,
                              "prob": low})

      upper_ps.append(local_probs_up)
      lower_ps.append(local_probs_low)
      local_probs_up = []
      local_probs_low = []
    return upper_ps, lower_ps

  start = 0
  finish = interval_length
  end = len(zscores)
  cont = True

  while cont:

    while start != finish:

      up, low = get_upper_lower_prob_from_score(score=zscores[start])

      local_probs_up.append({"idx":start,
                             "prob":up})

      local_probs_low.append({"idx":start,
                              "prob":low})
      start += 1
    upper_ps.append(local_probs_up)
    lower_ps.append(local_probs_low)
    local_probs_up = []
    local_probs_low =[]
    start = finish
    finish += interval_length

    if start >= len(zscores) or finish >= len(zscores):
      cont = False


  # pick up what is left
  while start != end:

    up, low = get_upper_lower_prob_from_score(score=zscores[start])

    local_probs_up.append({"idx":start,
                           "prob":up})

    local_probs_low.append({"idx":start,
                            "prob":low})
    start += 1

  if local_probs_low and local_probs_up:
    upper_ps.append(local_probs_up)
    lower_ps.append(local_probs_low)

  return upper_ps, lower_ps

def z_score_window_clusterer(windows, n_consecutive_windows, selector):
  '''
  classify each window according to the given states
  by computing a Z-score
  '''

  zscores = calculate_windows_zscore(windows=windows)
  upper_ps, lower_ps = windows_tails_p(zscores=zscores,
                                       interval_length=n_consecutive_windows)

  # we do the clustering for the windows
  for window in windows:
    selector(window=window,
             upper_ps=upper_ps,
             lower_ps=lower_ps)

  return windows


class ZScoreWindowCluster(object):

  def __init__(self, cutoff):
    self._cutoff = cutoff


  def __call__(self, window, upper_ps, lower_ps):

    window_upper_p, upper_vals = self.find_prob(widx=window.get_id(),
                                                probabilities=upper_ps)


    window_lower_p, lower_vals = self.find_prob(widx=window.get_id(),
                                                probabilities=lower_ps)

    # maximum upper p-value
    max_upper = np.amax(upper_vals)

    # maximum lower p-value
    max_lower = np.amax(lower_vals)

    # check if the p-value is such that
    # we reject the H0: NORMAL state and
    # Set the window to INSERT
    if max_upper < self._cutoff:
      window.set_state(state=WindowState.INSERT)
    elif max_lower < self._cutoff:

      # check if the p-value is such that
      # we reject the H0: NORMAL state and
      # Set the window to DELETE
      window.set_state(state=WindowState.DELETE)
    else:

      # the H0: NORMAL state was not rejected
      window.set_state(state=WindowState.NORMAL)

  def find_prob(self, widx, probabilities):

    for item in probabilities:
      for subitem in item:
        if subitem["idx"] == widx:
          return subitem["prob"], listify_dicts_property(list_dict_vals=item,
                                                         property_name="prob")

    raise Error("For window %s probability not found" % widx)

"""







