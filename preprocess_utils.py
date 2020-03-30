"""
Preprocessing utilities
"""

from pomegranate import*
import scipy.stats as st

from exceptions import Error
from helpers import windows_rd_statistics

VALID_DISTS = ['gaussian', 'uniform', 'poisson']


def fit_distribution(data, dist_name="gaussian", **kwargs):

    """
    Fits a distribution within the given dataset
    :param data:
    :param dist_name:
    :param kwargs:
    :return: appropriate distribution object
    """

    if dist_name not in VALID_DISTS:
        raise Error("Invalid distribution name. Name {0} not in {1}".format(dist_name, VALID_DISTS))

    if dist_name == 'gaussian':
        dist = NormalDistribution.from_samples(data)
        return dist
    elif dist_name == 'uniform':
        dist = UniformDistribution.from_samples(data)
        return dist
    elif dist_name == 'poisson':
        dist = PoissonDistribution.from_samples(data)
        return dist

def cluster(data, nclusters, method, **kwargs):

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
    clusterer.fit(data)
    return clusterer
  elif "zscore":
    return z_score_window_clusterer(windows=data,
                                    n_consecutive_windows=kwargs["n_consecutive_windows"],
                                    selector=kwargs["selector"])


  raise Error("Invalid clustering method: " + method )


def calculate_windows_zscore(windows):
  """
  Returns an array of Z-scores for every window
  in the given windows list

  Parameters
  ----------
  windows : list of Window instances
    DESCRIPTION.

  Returns
  -------
  list of floats representing the z-score for each window

  """

  mu = windows_rd_statistics(windows, statistic="mean")
  var = windows_rd_statistics(windows, statistic="var")

  z_scores = []

  for window in windows:
    zscore = (window.get_rd_stats(statistics="mean") - mu)/var
    z_scores.append(zscore)

  return z_scores

def windows_tails_p(zscores, l):
  """
  Calculate windows upper and lower tail probabilities arranged
  into window intervals of length l


  Parameters
  ----------
  zscores :  list
    DESCRIPTION. The score corresponding to each window
  l : int
    Length of window interval

  Returns
  -------
  upper_ps, lower_ps : list of upper probabilities,
                        list of lower probabilities

  """

  upper_ps = []
  lower_ps = []
  local_probs_up = []
  local_probs_low = []

  window_counter = 0
  for i, zscore in enumerate(zscores):

    if len(local_probs_up) < l:
      prob = st.norm.cdf(zscore)
      local_probs_up.append({"idx":i,
                             "prob":prob})

      local_probs_low.append({"idx":i,
                              "prob":1.0-prob})

    else:
      upper_ps.append(local_probs_up)
      lower_ps.append(local_probs_low)
      local_probs_up = []
      local_probs_low =[]

      prob = st.norm.cdf(zscore)
      prob = st.norm.cdf(zscore)
      local_probs_up.append({"idx":i,
                             "prob":prob})

      local_probs_low.append({"idx":i,
                              "prob":1.0-prob})

  return upper_ps, lower_ps



def z_score_window_clusterer(windows, n_consecutive_windows, selector):
  """
  classify each window according to the given states
  by computing a Z-score and a cut off quantity
  """

  zscores = calculate_windows_zscore(windows=windows)
  upper_ps, lower_ps = windows_tails_p(scores=zscores,
                                       l=n_consecutive_windows)

  # we do the clustering

  for window in windows:
    selector(window, upper_ps, lower_ps)

  return windows


class ZScoreWindowCluster(object):


  def __init__(self):
    pass

  def __call__(self, window, upper_ps, lower_ps):





