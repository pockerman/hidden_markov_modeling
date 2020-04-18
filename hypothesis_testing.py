"""
Hypothesis testing utilities
"""
from collections import defaultdict
import numpy as np
import scipy.stats as st

from preprocess_utils import compute_statistic
from helpers import WindowState

def zscore_statistic(data, null, **kwargs):


    statistic = compute_statistic(data=data,
                                  statistics = null.name)

    # compute the variance
    var = compute_statistic(data=data,
                                  statistics = "var")

    score = (statistic - null.value)/(np.sqrt(var))#/np.sqrt(len(data)))

    if "alternative" in kwargs:
      direction = kwargs["alternative"].direction

      if direction == ">" or direction == ">=":
        prob = 1.0 - st.norm.cdf(score)
        return prob, statistic
      elif direction == "<" or direction == "<=":
        prob = st.norm.cdf(score)
        return prob, statistic
      else:
        prob = 2.0*(1.0 - st.norm.cdf(score))
      return prob, statistic
    else:

      # assume two-sided by default
      prob = 2.0*(1.0 - st.norm.cdf(score))
      return prob, statistic

class Hypothesis(object):
  def __init__(self, name, direction, value):
    self._name = name
    self._direction = direction
    self._value = value
    self._description = name + direction +str(value)
    self._accepted=False

  @property
  def value(self):
    return self._value

  @value.setter
  def value(self, val):
    self._value = val

  @property
  def description(self):
    return self._description

  @property
  def direction(self):
    return self._direction

  @property
  def accepted(self):
    return self._accepted

  @accepted.setter
  def accepted(self, val):
    self._accepted = val

  @property
  def name(self):
    return self._name


class LessThan(Hypothesis):
  def __init__(self, parameter_name, value):
    super(LessThan, self).__init__(name=parameter_name,
                                   direction="<",
                                   value=value)


class LessOrEqualThan(Hypothesis):
  def __init__(self, parameter_name, value):
    super(LessOrEqualThan, self).__init__(name=parameter_name,
                                   direction="<=",
                                   value=value)

class GreaterThan(Hypothesis):
  def __init__(self, parameter_name, value):
    super(GreaterThan, self).__init__(name=parameter_name,
                                   direction=">",
                                   value=value)

class GreaterOrEqualThan(Hypothesis):
  def __init__(self, parameter_name, value):
    super(GreaterOrEqualThan, self).__init__(name=parameter_name,
                                   direction=">=",
                                   value=value)

class Equal(Hypothesis):
  def __init__(self, parameter_name, value):
    super(Equal, self).__init__(name=parameter_name,
                                   direction="=",
                                   value=value)

class NotEqual(Hypothesis):
  def __init__(self, parameter_name, value):
    super(NotEqual, self).__init__(name=parameter_name,
                                   direction="!=",
                                   value=value)


class HypothesisTest(object):

  def __init__(self, null, alternative,
               alpha, data, statistic_calculator):
    self._alpha = alpha
    self._null = null
    self._alternative = alternative
    self._data = data
    self._statistic_calculator = statistic_calculator
    self._p_value = None

  def test(self):
    self._p_value, statistic =self._statistic_calculator(data=self._data,
                                              null = self._null,
                                              alternative = self._alternative)

    print("Hypothesis test: ")
    if self._p_value < self._alpha :
      print("\tH0: " + self._null.description)
      print("\tvs")
      print("\tH1: " + self._alternative.description)
      print("\trejected H0 with a="+str(self._alpha) +
            " p-value="+str(self._p_value))
      print("\tstatistic computed: " +str(statistic))

      self._null.accepted=False
      self._alternative.accepted=True
    else:

      print("\tH0: " + self._null.description)
      print("\tvs")
      print("\tH1: " + self._alternative.description)
      print("\tcannot reject H0 with a="+str(self._alpha) +
            " p-value="+str(self._p_value))

      print("\tstatistic computed: " +str(statistic))
      self._null.accepted=True
      self._alternative.accepted=False


class SignificanceTestLabeler(object):

  def __init__(self, clusters, windows):
    self._clusters = clusters
    self._windows = windows

  def label(self, test_config):
    """
    Label the given clusters
    """

    cluster_data = defaultdict(list)

    labeled_clusters = {}

    if "NOT_NORMAL" in test_config["order"]:
       for cluster in self._clusters:

             print("Testing cluster: ", cluster.cidx)
             cluster_data[cluster.cidx] = cluster.get_data_from_windows(windows=self._windows)
             h0_normal = Equal(parameter_name=test_config["statistic_parameter"],
                                value=test_config["statistical_parameter_value"])

             ha_notnormal = NotEqual(parameter_name=test_config["statistic_parameter"],
                                                value=test_config["statistical_parameter_value"])
             hypothesis = HypothesisTest(null=h0_normal,
                                         alternative = ha_notnormal,
                                         alpha=test_config["significance"],
                                         data=cluster_data[cluster.cidx],
                                         statistic_calculator=zscore_statistic)
             hypothesis.test()
             if h0_normal.accepted:

                 print("Cluster %s is NORMAL" %cluster.cidx)
                 cluster.state = WindowState.NORMAL
                 if WindowState.NORMAL in labeled_clusters:
                    labeled_clusters[WindowState.NORMAL].merge(cluster)
                 else:
                    labeled_clusters[WindowState.NORMAL] = cluster

             else:

               cluster.state = WindowState.NOT_NORMAL
               if WindowState.NOT_NORMAL in labeled_clusters:
                 labeled_clusters[WindowState.NOT_NORMAL].merge(cluster)
               else:
                 labeled_clusters[WindowState.NOT_NORMAL] = cluster
    else:
      for cluster in self._clusters:

        print("Testing cluster: ", cluster.cidx)

        cluster_data[cluster.cidx] = cluster.get_data_from_windows(windows=self._windows)

        # is the cluster a DELETE or sth else:
        h0_delete = LessThan(parameter_name=test_config["statistic_parameter"],
                             value=test_config["statistical_parameter_value"])
        ha_delete = GreaterOrEqualThan(parameter_name=test_config["statistic_parameter"],
                                       value=test_config["statistical_parameter_value"])

        hypothesis_delete = HypothesisTest(null=h0_delete,
                                           alternative = ha_delete,
                                           alpha=test_config["significance"],
                                           data=cluster_data[cluster.cidx],
                                           statistic_calculator=zscore_statistic)

        hypothesis_delete.test()

        if h0_delete.accepted:
        # this cluster is a delete cluster
          print("Cluster %s is DELETE" %cluster.cidx)
          cluster.state = WindowState.DELETE

          if WindowState.DELETE in labeled_clusters:
                labeled_clusters[WindowState.DELETE].merge(cluster)
          else:
                labeled_clusters[WindowState.DELETE] = cluster

        else:

          h0_normal = Equal(parameter_name=test_config["statistic_parameter"],
                             value=test_config["statistical_parameter_value"])
          ha_normal = GreaterThan(parameter_name=test_config["statistic_parameter"],
                             value=test_config["statistical_parameter_value"])

          hypothesis_normal = HypothesisTest(null=h0_normal,
                                             alternative = ha_normal,
                                             alpha=test_config["significance"],
                                             data=cluster_data[cluster.cidx],
                                             statistic_calculator=zscore_statistic)

          hypothesis_normal.test()

          if h0_normal.accepted:
            print("Cluster %s is NORMAL"%cluster.cidx)
            cluster.state = WindowState.NORMAL

            if WindowState.NORMAL in labeled_clusters:
              labeled_clusters[WindowState.NORMAL].merge(cluster)
            else:
              labeled_clusters[WindowState.NORMAL] = cluster
          else:
            print("Cluster %s is INSERT"%cluster.cidx)
            cluster.state = WindowState.INSERT

            if WindowState.INSERT in labeled_clusters:
              labeled_clusters[WindowState.INSERT].merge(cluster)
            else:
              labeled_clusters[WindowState.INSERT] = cluster

    return labeled_clusters







