import json
import numpy as np
import logging
import json
from enum import Enum
from scipy import stats

from collections import namedtuple
from exceptions import FullWindowException
from exceptions import Error

DUMMY_ID = -1
INFO="INFO:"
WARNING="WARNING:"

class HMMCallback(object):

  def __init__(self, callback):
    self._callback = callback
    self.model = None

  def on_epoch_end(self, logs):
    pass
    #self._callback(logs)

  def on_training_begin(self):
    pass

  def on_training_end(self, logs):
    pass

def print_logs_callback(logs):
  print(logs)

def read_configuration_file(config_file):
    """
    Read the json configuration file and
    return a map with the configuration
    entries
    """
    with open(config_file) as json_file:
        configuration = json.load(json_file)
        return configuration


def save_windows(windows, configuration, win_interval_length):

  if configuration["save_windows"]:
    import json
    with open(configuration["windows_filename"]+
                  "_"+str(win_interval_length)+".json", 'w') as jsonfile:
      json_str = windows_to_json(windows)
      json.dump(json_str, jsonfile)

def save_hmm(hmm_model, configuration, win_interval_length):

  if configuration["HMM"]["save_model"]:
    json_str = hmm_model.to_json()
    import json
    with open(configuration["HMM"]["save_hmm_filename"]+
              "_"+str(win_interval_length)+".json", 'w') as jsonfile:
      json.dump(json_str, jsonfile)

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

def flat_windows_from_state(windows, configuration, as_on_seq):

  """
  Returns a flattened list of windows by their
  prop property observations
  """

  win = []

  for window in windows:

    if not isinstance(window, Window):
      raise Error("The given window is not an insatnce of Window")

    if window.get_state() == WindowState.NORMAL:
      value = [configuration["normal_rd"]] if as_on_seq else configuration["normal_rd"]
      win.append(value)
    elif window.get_state() == WindowState.DELETE:
      value = [configuration["delete_rd"]] if as_on_seq else configuration["delete_rd"]
      win.append(value)
    elif window.get_state() == WindowState.INSERT:
      value = [configuration["insert_rd"]] if as_on_seq else configuration["insert_rd"]
      win.append(value)

  return win


def flat_windows_rd_from_indexes(indexes, windows):
  rd_observations = []

  if indexes is None:
    # do all the windows
    for window in windows:
      rd_observations.extend(window.get_rd_observations())
  else:

    for widx in indexes:
      rd_observations.extend(windows[widx].get_rd_observations())
  return rd_observations

"""
def windows_rd_statistics(windows, statistic="all"):

  rd_observations = []

  for window in windows:
    rd_observations.extend(window.get_rd_observations())

  return compute_statistic(data=rd_observations, statistsics=statistic)
"""


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


def add_window_observation(window, windows,
                           observation, windowcapacity):
    """
    Add a new observation to the given window. If the
    window has reached its capacity a new window
    is created and then the observation is appened
    :param window: The window instance to add the observation
    :param windows: The list of windows where the window is cached
    :param observation: The observation to add in the window
    :param windowcapacity: The maximum window capacity
    :return: instance of Window class
    """

    if window.has_capacity():
        window.add(observation=observation)
    else:
        windows.append(window)
        window = Window(idx=window.get_id()+1,
                        capacity=windowcapacity)
        window.add(observation=observation)

    return window


class WindowState(Enum):
  DELETE = 0
  NORMAL = 1
  INSERT = 2
  TUF = 3
  NOT_NORMAL = 4
  INVALID = 5


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

        # holds tuples of observations:
        # the first tuple is the wga_treated
        # and the second is the non wga_treated
        self._observations = []

        # the total number of insertions/deletions
        self._net_indels = 0

        # the state of the window
        self._state = WindowState.INVALID

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
                            "median", "min", "max",
                            "mode"]

        if statistics not in valid_statistics:
            raise Error("Invalid statistsics: '{0}'"
                        " not in {1}".format(statistics, valid_statistics))

        # accumulate RD as an array and use numpy
        rd_data = [item.read_depth for item in self._observations]
        from preprocess_utils import compute_statistic
        return compute_statistic(data=rd_data,statistics=statistics)

    def get_rd_observations(self):
        """
        Returns a list with read depth observations
        contained in the window
        :return:
        """
        rd_data = [item.read_depth for item in self._observations]
        return rd_data

    def get_bases(self):

      bases = []

      for observation in self._observations:
        bases.append(observation.base[0])

      return bases

    def get_sequence(self):

      seq = ""
      for observaion in self._observations:
        seq += observaion.base[0]
      return seq.strip()


    def get_gc_percent(self):
      gc_count = 0

      for observation in self._observations:
         if observation.base[0].upper() == "C" \
                    or observation.base[0].upper() == "G":
                gc_count += 1

      if gc_count == 0:
        return 0

      return gc_count / len(self._observations)


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
      """
      Set the state of the window

      Parameters
      ----------
      state : WindowState
        Enumeration describing the window state

      Returns
      -------
      None.

      """
      self._state = state

    def get_state(self):
      """
      the state of the window

      Returns
      -------
      WindowState

      """
      return self._state

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
                  "state":self._state.name,
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


class MixedWindowView(object):
  """
  A class that holds two instances of windows
  """

  def __init__(self, wga_w, n_wga_w):
    self._windows={"wga_w": wga_w,
                   "n_wga_w": n_wga_w}

    # the state of the window
    self._state = WindowState.INVALID


  @property
  def state(self):
    return self._state

  @state.setter
  def state(self, value):
    self._state = value

  def get_rd_counts(self, name):

     if name == "both":
          return (self._windows["wga_w"].get_rd_observations(), self._windows["n_wga_w"].get_rd_observations())
     elif name == 'wga_w':
          return self._windows["wga_w"].get_rd_observations()
     elif name == 'n_wga_w':
          return self._windows["n_wga_w"].get_rd_observations()
     else:
          raise Error("Name {0} is invalid ".format(name))


  def get_rd_stats(self, statistics="all", name="both"):
        """
        Returns a statistical summary as a dictionary
        of the read depth variable in the window
        :param statistics:
        :return:
        """

        if name == "both":
          return (self._windows["wga_w"].get_rd_stats(statistics=statistics), self._windows["n_wga_w"].get_rd_stats(statistics=statistics))
        elif name == 'wga_w':
          return self._windows["wga_w"].get_rd_stats(statistics=statistics)
        elif name == 'n_wga_w':
          return self._windows["n_wga_w"].get_rd_stats(statistics=statistics)
        else:
          raise Error("Name {0} is invalid ".format(name))

  def get_bases(self, windowtype="both"):

    if windowtype == "both":
      wga_w = self._windows["wga_w"]
      n_wga_w = self._windows["n_wga_w"]

      pairs = zip(wga_w.get_bases(), n_wga_w.get_bases())
      return pairs

  def get_sequence(self, windowtype="both"):

    if windowtype == "both":
      return self._windows["wga_w"].get_sequence(),\
        self._windows["n_wga_w"].get_sequence()

  def get_gc_percent(self, windowtype="both"):
    if windowtype == "both":
      return self._windows["wga_w"].get_gc_percent(),\
        self._windows["n_wga_w"].get_gc_percent()










