import json
import numpy as np
import logging
import json
from enum import Enum

from collections import namedtuple
from exceptions import FullWindowException
from exceptions import Error

DUMMY_ID = -1

def read_configuration_file(config_file):
    """
    Read the json configuration file and
    return a map with the configuration
    entries
    """
    with open(config_file) as json_file:
        configuration = json.load(json_file)
        return configuration


def listify_dicts_property(list_dict_vals, property_name):
  """
  given a list of dictionaries return a list with the
  values of the property with name property_name

  Parameters
  ----------
  list_dict_vals : list of dictioanries

  property_name : str
    The property name to look for

  Returns
  -------
  result : list


  """

  result = []

  for item in list_dict_vals:
    if item.get(property_name):
      result.append(item[property_name])

  return result


def set_up_logger(configuration):
    # set up the logger
    logger_file = configuration.get("logger_file", None)

    # always force logging
    if logger_file is None:
        logger_file = "tuf.log"

    logging_level = configuration.get("logger_level", None)

    if logging_level is None or logging_level == "INFO":
        logging_level = logging.INFO
    elif logging_level == "DEBUG":
        logging_level = logging.DEBUG
    elif logging_level == "ERROR":
        logging_level = logging.ERROR

    logging.basicConfig(filename=logger_file, level=logging_level)

def flat_windows(windows, prop="RD"):

  """
  Returns a flattened list of windows by their
  prop property observations
  """

  win = []

  for window in windows:

    if not isinstance(window, Window):
      raise Error("The given window is not an insatnce of Window")

    win.append(window.get_rd_observations())

  return win

def windows_rd_statistics(windows, statistic="all"):

  rd_observations = []

  for window in windows:
    rd_observations.extend(window.get_rd_observations())

  if statistic == "mean":
    return np.mean(rd_observations)
  elif statistic == "var":
    return np.var(rd_observations)
  elif statistic == "median":
    return np.median(rd_observations)
  elif statistic == "all":
    mu = np.mean(rd_observations)
    var = np.var(rd_observations)
    median = np.median(rd_observations)
    return mu, var, median
  else:
    raise Error("Unknown statistic: %s" % statistic)


def windows_to_json(windows):
  """
  Generate a json formatted string
  representing the windows
  Parameters
  ----------
  windows : list
    DESCRIPTION. A list of Window instances

  Returns
  -------
  string

  """

  str_map = {}

  for window in windows:
    str_map[window.get_id()] = window.to_json()

  return json.dumps(str_map)


def windows_from_json(jsonmap):
  """


  Parameters
  ----------
  jsonmap : TYPE
    DESCRIPTION.

  Returns
  -------
  list of Window instances

  """

  windows = []

  for idx in jsonmap.keys():

    window_str = jsonmap[idx]
    window_data = json.loads(window_str)
    window = Window(idx=window_data["id"],
                    capacity=window_data["capacity"])

    for obs in window_data["observations"]:
      obsdata = json.loads(obs)
      observaion =observation_from_json(jsonmap=obsdata)
      window.add(observation=observaion)

    windows.append(window)
  return windows


Observation = namedtuple("Observation", ["position", "read_depth", "base"])

def observation_to_json(observation):
  """
  Returns a json formatted string for the
  given observation

  Parameters
  ----------
  observation : Observation
    DESCRIPTION.

  Returns
  -------
  json_str : string

  """

  json_str = {"position":observation.position,
              "read_depth": observation.read_depth,
              "base": observation.base}

  json_str = json.dumps(json_str)
  return json_str


def observation_from_json(jsonmap):

  observation = Observation(position=jsonmap["position"],
                            read_depth=jsonmap["read_depth"],
                            base=jsonmap["base"])
  return observation


class WindowState(Enum):
  DELETE = 0
  NORMAL = 1
  INSERT = 2
  TUF = 3


class WindowIterator(object):

    """
    Helper class to allow iteration over the window
    elements
    """

    def __init__(self, data):

        # data over which to iterate
        self._data = data

        # current index
        self._counter = 0

    def __next__(self):
        if self._counter < len(self._data):
            tmp = self._data[self._counter]
            self._counter += 1
            return tmp

        raise StopIteration


class Window(object):
    """
    Class representing a window for arranging of the data
    """
    def __init__(self, idx, capacity):

        # the id of the window
        self._id = idx

        # maximum capacity of the window
        self._capacity = capacity

        # holds the observations i.e. the base
        # strings observed into the window
        self._observations = []

        # the total number of insertions/deletions
        self._net_indels = 0

        # the state of the window
        self._state = None

    def add(self, observation):

        if len(self._observations) >= self._capacity:
            raise FullWindowException(self._capacity)

        # all the observations accumulated in the window
        self._observations.append(observation)

    def capacity(self):
        """
        Returns the capacity of the window i.e. the maximum
        number of observations the window can accumulate
        :return:
        """
        return self._capacity

    def has_capacity(self):
        """
        Returns True if the caapacity of the
        window has not been exhausted
        :return: boolean
        """
        return (self.capacity() - len(self._observations)) != 0

    def get_range(self, start, end):
        """
        Returns a range of the observations stored in
        the window
        :param start:
        :param end:
        :return:
        """
        return self._observations[start:end]

    def get_rd_stats(self, statistics="all"):
        """
        Returns a statistical summary as a dictionary
        of the read depth variable in the window
        :param statistics:
        :return:
        """
        valid_statistics = ["all",  "mean", "var",
                            "median", "min", "max"]

        if statistics not in valid_statistics:
            raise Error("Invalid statistsics: '{0}'"
                        " not in {1}".format(statistics, valid_statistics))

        # accumulate RD as an array and use numpy
        rd_data = [item.read_depth for item in self._observations]

        if statistics == "mean":
          return np.mean(rd_data)
        elif statistics == "var":
          return np.var(rd_data)
        elif statistics == "median":
          return np.median(rd_data)
        elif statistics == "min":
          return np.amin(rd_data)
        elif statistics == "max":
          return np.amax(rd_data)
        elif statistics == "all":

          mean = np.mean(rd_data)
          var = np.var(rd_data)
          median = np.median(rd_data)
          min = np.amin(rd_data)
          max = np.amax(rd_data)
          return {"mean": mean, "var": var,
                  "median": median, "min": min,
                  "max": max}

    def get_rd_observations(self):
        """
        Returns a list with read depth observations
        contained in the window
        :return:
        """
        rd_data = [item.read_depth for item in self._observations]
        return rd_data

    def get_gc_count(self):
        """
        Returns the GC count for the window
        :return:
        """
        gc_count = 0
        at_count = 0

        for observation in self._observations:

            if observation.base.upper() == "C" \
                    or observation.base.upper() == "G":
                gc_count += 1
            elif observation.base.upper() != "N":
                at_count += 1

        if gc_count == 0:
            return 0

        return gc_count / (gc_count + at_count)

    def set_net_indels(self, net_indels):
        """
        Set the net insertion/deletions for the window
        :param net_indels:
        :return:
        """
        self._net_indels = net_indels

    def has_gaps(self):
        """
        Returns true if the window contains gaps. This is done
        by looping over the the window observations and check if the
        observed positions are contiguous without gaps
        :return:
        """
        previous = int(self._observations[0].position)
        for item in range(1, len(self)):

            pos = int (self._observations[item].position)

            if pos != previous + 1:
                # we have a gap
                return True
            else:
                previous = pos
        return False


    def get_id(self):
      """
      Returns the zero based id of the window

      Returns
      -------
      self._id

      """
      return self._id

    def set_state(self, state):
      self._state = state

    def insert_at(self, pos, data):
        """
        Insert the data at the specified position
        for the current window
        :param pos:
        :param data:
        :return:
        """
        self._observations.insert(pos, data)

    def to_json(self):
      """
      Returns a json formatted string represneting
      the window
      """
      json_str = {"id":self._id,
                  "capacity":self._capacity,}

      observations = []
      for obs in self._observations:
        obs_json = observation_to_json(observation=obs)
        observations.append(obs_json)

      json_str["observations"] = observations
      return json.dumps(json_str)

    def __len__(self):
        return len(self._observations)

    def __iter__(self):
        """
        Produce an iterator to iterate over the accumulated
        window observations
        :return:
        """
        return WindowIterator(data=self._observations)

    def __getitem__(self, item):
        """
        Returns the item-th observation
        :param item: int the index of the observation to access
        :return: Observation instance
        """
        return self._observations[item]

    def __setitem__(self, o, value):
        """
        Set the o-th observation to value
        :param o:
        :param value:
        :return:
        """
        self._observations[o] = value


