"""
Hypothesis testing utilities
"""
from collections import defaultdict
import numpy as np
import scipy.stats as st

from preprocess_utils import compute_statistic
from helpers import WindowState
from helpers import INFO

def zscore_statistic(data, null, **kwargs):


    statistic = compute_statistic(data=data,
                                  statistics = null.name)

    # compute the variance
    var = compute_statistic(data=data, statistics = "var")

    score = (statistic - null.value)/(np.sqrt(var))

    if "alternative" in kwargs:
      direction = kwargs["alternative"].direction

      if direction == ">" or direction == ">=":
        prob = 1.0 - st.norm.cdf(score)
        return prob, statistic
      elif direction == "<" or direction == "<=":
        prob = st.norm.cdf(score)
        return prob, statistic
      else:
        prob = 2.0*(1.0 - st.norm.cdf(np.abs(score)))
        return prob, statistic
    else:

      # assume two-sided by default
      prob = 2.0*(1.0 - st.norm.cdf(np.fabs(score)))
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
    self._p_value, statistic = \
      self._statistic_calculator(data=self._data,
                                 null = self._null,
                                 alternative = self._alternative)

    print("{0} Hypothesis test: ".format(INFO))
    if self._p_value < self._alpha :
      print("\tH0: " + self._null.description +
            "vs" + "H1: " + self._alternative.description)

      print("\trejected H0 with a="+str(self._alpha) +
            " p-value="+str(self._p_value))
      print("\tstatistic computed: " +str(statistic))

      self._null.accepted=False
      self._alternative.accepted=True
    else:

      print("\tH0: " + self._null.description +
            "vs" + "H1: " + self._alternative.description)
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

    print("{0} Labeling clusters...".format(INFO))
    cluster_data = defaultdict(list)


    for cluster in self._clusters:

        print("{0} Testing cluster {1}".format(INFO, cluster.cidx))
        print("{0} Testing FULL DELETE".format(INFO))

        windows = [window.get_window(wtype="wga_w")
                     for window in self._windows]

        cluster_data = cluster.get_data_from_windows(windows=windows)

        # get the cluster statistics
        cluster_stats = cluster.get_statistics(windows=self._windows,
                                               statistic="all",
                                               wtype="wga_w")

        h0 = \
          Equal(parameter_name=test_config["statistic_parameter"],
                value=0.0)

        ha = \
          GreaterThan(parameter_name=test_config["statistic_parameter"],
                      value=0.0)

        test = HypothesisTest(null=h0,
                              alternative=ha,
                              alpha=test_config["significance"],
                              data=cluster_data,
                              statistic_calculator=zscore_statistic)

        test.test()

        if h0.accepted:
          print("{0} Cluster {1} is  labeled as DELETE".format(INFO,cluster.cidx))
          cluster.state = WindowState.DELETE
        else:

          print("{0} Testing ONE COPY DELETE".format(INFO))
          h0 = \
            Equal(parameter_name=test_config["statistic_parameter"],
                  value=10.0)

          ha = \
              NotEqual(parameter_name=test_config["statistic_parameter"],
                          value=10.0)

          test = HypothesisTest(null=h0,
                                alternative=ha,
                                alpha=test_config["significance"],
                                data=cluster_data[cluster.cidx],
                                statistic_calculator=zscore_statistic)

          test.test()
          if h0.accepted:
            print("{0} Cluster {1} is  labeled as ONE COPY DELETE".format(INFO, cluster.cidx))
            cluster.state = WindowState.ONE_COPY_DELETE
          else:
            print("{0} Testing NORMAL".format(INFO))

            h0 = \
              Equal(parameter_name=test_config["statistic_parameter"],
                  value=20.0)

            ha = \
                NotEqual(parameter_name=test_config["statistic_parameter"],
                          value=20.0)

            test = HypothesisTest(null=h0,
                                  alternative=ha,
                                  alpha=test_config["significance"],
                                  data=cluster_data[cluster.cidx],
                                  statistic_calculator=zscore_statistic)

            test.test()
            if h0.accepted:
              print("{0} Cluster {1} is  labeled as NORMAL".format(INFO, cluster.cidx))
              cluster.state = WindowState.NORMAL
            else:

              print("{0} Testing for DUPLICATION".format(INFO))

              h0 = \
                GreaterThan(parameter_name=test_config["statistic_parameter"],
                            value=20.0)

              ha = \
                LessThan(parameter_name=test_config["statistic_parameter"],
                          value=20.0)

              test = HypothesisTest(null=h0,
                                      alternative=ha,
                                      alpha=test_config["significance"],
                                      data=cluster_data[cluster.cidx],
                                      statistic_calculator=zscore_statistic)

              test.test()
              if h0.accepted:
                print("{0} Cluster {1} is  labeled as DUPLICATION".format(INFO, cluster.cidx))
                cluster.state = WindowState.INSERT
              else:
                # if we reach here then this means that statistically
                # mu is not 10 is not 20 is not 0 and it is
                # less than 20 we classify this cluster as TUF
                print("{0} Cluster {1} is  labeled as TUF".format(INFO, cluster.cidx))
                cluster.state = WindowState.TUF


    print("{0} Done...".format(INFO))
    return self._clusters







