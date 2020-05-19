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
ERROR="ERROR:"
DEBUG="DEBUG:"


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


class Observation(object):
  def __init__(self, position, read_depth, base):
    self._position = int(position)
    self._read_depth = int(read_depth)

    if isinstance(base, str):
      self._base = [base]
    elif isinstance(base, list):
      self._base = base
    else:
      raise Error("Unknown type for base "
                  "in observation. Type {1} "
                  "not in ['str', 'list']".format(type(base)))

  @property
  def position(self):
    return self._position

  @position.setter
  def position(self, value):
    self._position = value

  @property
  def read_depth(self):
    return self._read_depth

  @read_depth.setter
  def read_depth(self, value):
    self._read_depth = value

  @property
  def base(self):
    return self._base

  @base.setter
  def base(self, value):
    return self._base


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
        window = Window(idx=window.idx + 1,
                        capacity=windowcapacity)
        window.add(observation=observation)

    return window


class WindowType(Enum):
  WGA = 0
  NO_WGA = 1
  BOTH = 2
  N_WIN = 3

  @staticmethod
  def from_string(string):
    if string.upper() == "WGA":
      return WindowType.WGA
    elif string.upper() == "NO_WGA":
      return WindowType.NO_WGA
    elif string.upper() == "BOTH":
      return WindowType.BOTH
    elif string.upper() == 'N_WIN':
      return WindowType.N_WIN

    raise Error("Invalid WindowType. "
                "Type {0} not in {1}".format(string,
                                             ["WGA",
                                              "NO_WGA",
                                              "BOTH",
                                              "N_WIN"]))

class WindowState(Enum):
  DELETE = 0
  ONE_COPY_DELETE = 1
  NORMAL = 2
  INSERT = 3
  TUF = 4
  NOT_NORMAL = 5
  OTHER = 6
  INVALID = 7

  @staticmethod
  def from_string(string):
    if string.upper() == "DELETE":
      return WindowState.DELETE
    elif string.upper() == "ONE_COPY_DELETE":
      return WindowState.ONE_COPY_DELETE
    elif string.upper() == "NORMAL":
      return WindowState.NORMAL
    elif string.upper() == "INSERT":
      return WindowState.INSERT
    elif string.upper() == "TUF":
      return WindowState.TUF
    elif string.upper() == "OTHER":
      return WindowState.OTHER

    raise Error("Invalid WindowState. "
                "Type {0} not in {1}".format(string,
                                             ["DELETE",
                                              "ONE_COPY_DELETE",
                                              "NORMAL",
                                              "INSERT", "TUF", "OTHER"]))



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

    N_WINDOW_MARKER=-999

    @staticmethod
    def set_window_marker(marker):
        Window.N_WINDOW_MARKER= marker

    def __init__(self, idx, capacity):

        # the id of the window
        self._id = idx

        # maximum capacity of the window
        self._capacity = capacity

        # holds tuples of observations:
        # the first tuple is the wga_treated
        # and the second is the non wga_treated
        self._observations = []

        # the state of the window
        self._state = WindowState.INVALID

    @property
    def idx(self):
      return self._id

    @property
    def state(self):
      return self._state

    @state.setter
    def state(self, value):
      self._state = value

    @property
    def capacity(self):
        """
        Returns the capacity of the window i.e. the maximum
        number of observations the window can accumulate
        :return:
        """
        return self._capacity

    def get_start_end_pos(self):
      return (self._observations[0].position,
              self._observations[len(self._observations)-1].position)

    def has_capacity(self):
        """
        Returns True if the caapacity of the
        window has not been exhausted
        :return: boolean
        """
        return (self.capacity - len(self._observations)) != 0

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
                        " not in {1}".format(statistics,
                                             valid_statistics))

        # accumulate RD as an array and use numpy
        rd_data = [item.read_depth for item in self._observations]
        from preprocess_utils import compute_statistic
        return compute_statistic(data=rd_data,
                                 statistics=statistics)


    def add(self, observation):

        if len(self._observations) >= self._capacity:
            raise FullWindowException(self._capacity)

        # all the observations accumulated in the window
        self._observations.append(observation)

    def set_window_rd_mark(self, mark):
      try:
        Window.set_window_marker(marker=mark)
        for obs in self._observations:
          try:
            obs.read_depth=mark
          except Exception as e:
            print("{0} At observation {1} {2}".format(ERROR, obs, str(e)))
            raise
      except Exception as e:
        print("{0} For window {1} an excpetion {2}".format(ERROR, self.idx, str(e)))
        raise

    def get_rd_observations(self):
        """
        Returns a list with read depth observations
        contained in the window
        :return:
        """
        rd_data = [item.read_depth
                   for item in self._observations]
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

    def has_base(self, string):
      for observaion in self._observations:
        if observaion.base[0].upper() == string.upper():
          return True
      return False

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
    self._windows={WindowType.WGA: wga_w,
                   WindowType.NO_WGA: n_wga_w}

    # the state of the window
    self._state = WindowState.INVALID


  @property
  def state(self):
    return self._state

  @state.setter
  def state(self, value):
    self._state = value


  def is_n_window(self):

    if self._windows[WindowType.WGA].state == WindowType.N_WIN or\
      self._windows[WindowType.NO_WGA].state == WindowType.N_WIN:
        return True

    return False

  def get_window(self, wtype):
    return self._windows[wtype]

  def get_rd_counts(self, name):

     if name == WindowType.BOTH:
          return (self._windows[WindowType.WGA].get_rd_observations(),
                  self._windows[WindowType.NO_WGA].get_rd_observations())
     elif name == WindowType.WGA:
          return self._windows[WindowType.WGA].get_rd_observations()
     elif name == WindowType.NO_WGA:
          return self._windows[WindowType.NO_WGA].get_rd_observations()

     raise Error("Name {0} is invalid ".format(name))


  def get_rd_stats(self, statistics="all", name=WindowType.BOTH):
        """
        Returns a statistical summary as a dictionary
        of the read depth variable in the window
        :param statistics:
        :return:
        """
        if self.is_n_window():
           if name == WindowType.BOTH:
             return (Window.N_WINDOW_MARKER, Window.N_WINDOW_MARKER)
           else:
            return Window.N_WINDOW_MARKER

        if name == WindowType.BOTH:
          return (self._windows[WindowType.WGA].get_rd_stats(statistics=statistics),
                  self._windows[WindowType.NO_WGA].get_rd_stats(statistics=statistics))
        elif name == WindowType.WGA:
          return self._windows[WindowType.WGA].get_rd_stats(statistics=statistics)
        elif name == WindowType.NO_WGA:
          return self._windows[WindowType.NO_WGA].get_rd_stats(statistics=statistics)

        raise Error("Windowtype {0}"
                    " not in {1}".format(name, [WindowType.BOTH.name,
                                                WindowType.WGA.name,
                                                WindowType.NO_WGA.name]))

  def get_bases(self, windowtype="both"):

    if windowtype == "both":
      wga_w = self._windows[WindowType.WGA]
      n_wga_w = self._windows[WindowType.NO_WGA]

      pairs = zip(wga_w.get_bases(), n_wga_w.get_bases())
      return pairs

  def get_sequence(self, windowtype="both"):

    if windowtype == "both":
      return (self._windows[WindowType.WGA].get_sequence(),
              self._windows[WindowType.NO_WGA].get_sequence())
    elif windowtype == WindowType.WGA:
      return self._windows[WindowType.WGA].get_sequence()
    elif windowtype == WindowType.NO_WGA:
      return self._windows[WindowType.NO_WGA].get_sequence()


    raise Error("Windowtype {0}"
                " not in {1}".format(windowtype, ["both",
                                                  WindowType.WGA.name,
                                                  WindowType.NO_WGA.name]))

  def get_gc_percent(self, windowtype="both"):
    if windowtype == "both":
      return self._windows[WindowType.WGA].get_gc_percent(),\
        self._windows[WindowType.NO_WGA].get_gc_percent()
    elif windowtype == WindowType.WGA:
      return self._windows[WindowType.WGA].get_gc_percent()
    elif windowtype == WindowType.NO_WGA:
      return self._windows[WindowType.NO_WGA].get_gc_percent()

    raise Error("Windowtype {0}"
                " not in {1}".format(windowtype, ["both",
                                                  WindowType.WGA.name,
                                                  WindowType.NO_WGA.name]))













