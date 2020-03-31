"""
Preprocessing utilities
"""

from pomegranate import*
import scipy.stats as st
import numpy as np

from exceptions import Error
from helpers import windows_rd_statistics
from helpers import listify_dicts_property
from helpers import WindowState

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

def windows_tails_p(zscores, interval_length):
  """
  Calculate windows upper and lower tail probabilities arranged
  into window intervals of length l

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

  """

  local_probs_up = []
  local_probs_low = []
  upper_ps = []
  lower_ps = []

  if interval_length == 1:
    for i, zscore in enumerate(zscores):

      prob = st.norm.cdf(zscore)
      local_probs_up.append({"idx":i,
                             "prob":prob})

      local_probs_low.append({"idx":i,
                              "prob":1.0-prob})

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
      prob = st.norm.cdf(zscores[start])
      local_probs_up.append({"idx":start,
                             "prob":prob})

      local_probs_low.append({"idx":start,
                              "prob":1.0-prob})
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

    prob = st.norm.cdf(zscores[start])
    local_probs_up.append({"idx":start,
                             "prob":prob})

    local_probs_low.append({"idx":start,
                              "prob":1.0-prob})
    start += 1

  if local_probs_low and local_probs_up:
    upper_ps.append(local_probs_up)
    lower_ps.append(local_probs_low)

  #we also need to last ones
  return upper_ps, lower_ps

def z_score_window_clusterer(windows, n_consecutive_windows, selector):
  """
  classify each window according to the given states
  by computing a Z-score and a cut off quantity
  """

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

    max_upper = np.amax(upper_vals)
    max_lower = np.amax(lower_vals)

    if max_upper < self._cutoff:
      window.set_state(state=WindowState.INSERT)
    elif max_lower < self._cutoff:
      window.set_state(state=WindowState.DELETE)
    else:
      window.set_state(state=WindowState.NORMAL)

  def find_prob(self, widx, probabilities):

    for item in probabilities:
      for subitem in item:
        if subitem["idx"] == widx:
          return subitem["prob"], listify_dicts_property(list_dict_vals=item,
                                                    property_name="prob")

    raise Error("For window %s probability not found" % widx)









