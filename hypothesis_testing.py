"""
Hypothesis testing utilities
"""

import numpy as np
import scipy.stats as st

from preprocess_utils import compute_statistic

def zscore_statistic(data, null, **kwargs):


    statistic = compute_statistic(data=data,
                                  statistics = null.name)

    # compute the variance
    var = compute_statistic(data=data,
                                  statistics = "var")

    score = (statistic - null.value)/(np.sqrt(var)/np.sqrt(len(data)))

    if "alternative" in kwargs:
      direction = kwargs["alternative"].direction

      if direction == ">" or direction == ">=":
        prob = 1.0 - st.norm.cdf(score)
        return prob, statistic
      elif direction == ">" or direction == ">=":
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






