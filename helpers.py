import json
import logging
import time
from enum import Enum
from functools import wraps

from exceptions import Error

DUMMY_ID = -1
INFO="INFO:"
WARNING="WARNING:"
ERROR="ERROR:"
DEBUG="DEBUG:"

def partition_range(start, end, npieces):

  if npieces == 0:
    raise Error("Zero number of partitions")

  load = (end - start )// npieces
  chunks = []

  start_p = start

  # end is inclusive
  end_p = start_p + load -1
  for p in range(npieces-1):

    chunks.append((start_p, end_p ))
    start_p = end_p + 1
    end_p += load

  if end_p != end:
    end_p += end - end_p

  #remain = (end - start ) % npieces
  #if remain != 0:
  #  end_p += remain

  chunks.append((start_p, end_p ))

  return chunks

def read_configuration_file(config_file):
    """
    Read the json configuration file and
    return a map with the config entries
    """
    with open(config_file) as json_file:
        configuration = json.load(json_file)
        return configuration

def timefn(fn):
  @wraps(fn)
  def measure(*args, **kwargs):
    time_start = time.perf_counter()
    result = fn(*args, **kwargs)
    time_end = time.perf_counter()
    print("{0} Done. Execution time"
          " {1} secs".format(INFO, time_end - time_start))
    return result
  return measure


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
  INSERTION = 3
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
    elif string.upper() == "INSERTION":
      return WindowState.INSERTION
    elif string.upper() == "TUF":
      return WindowState.TUF
    elif string.upper() == "OTHER":
      return WindowState.OTHER

    raise Error("Invalid WindowState. "
                "Type '{0}' not in {1}".format(string,
                                             ["DELETE",
                                              "ONE_COPY_DELETE",
                                              "NORMAL",
                                              "INSERTION", "TUF", "OTHER"]))


class Window(object):
    """
    Class representing a window for arranging of the data
    """

    __slots__ = ('_id', '_capacity', '_samdata', '_state' )

    N_WINDOW_MARKER=-999


    @staticmethod
    def set_window_marker(marker):
        Window.N_WINDOW_MARKER= marker

    def __init__(self, idx, capacity, samdata):

        # the id of the window
        self._id = idx

        # maximum capacity of the window
        self._capacity = capacity

        # the data collected from the SAM file
        self._samdata = samdata

        # the state of the window
        self._state = WindowState.INVALID

    @property
    def idx(self):
      return self._id

    @property
    def start(self):
      return self._samdata["start"]

    @property
    def end(self):
      return self._samdata["end"]

    @property
    def state(self):
      return self._state

    @state.setter
    def state(self, value):
      self._state = value

    @property
    def capacity(self):
        return self._capacity

    @property
    def start_end_pos(self):
      return self.sam_property("start"), self.sam_property("end")

    def sam_property_names(self):
      return self._samdata.keys()

    def sam_property(self, name):
      return self._samdata[name]

    def get_rd_statistic(self, statistic):

      if statistic== "mean":
        return self._samdata["qmean"]

      elif statistic == "median":
        self._samdata["qmedian"]

      raise Error("Statistic '{0}' is not currently "
                  "computed for Window".format(statistic))

    def get_gc_percent(self):
     return self._samdata["gcr"]

    def set_window_rd_mark(self, mark):
        Window.set_window_marker(marker=mark)

    def has_gaps(self):
        """
        Returns true if the window contains gaps.
        """
        return self._samdata["gapAlert"]



class MixedWindowView(object):
  """
  A class that holds two instances of windows
  """

  __slots__ = ('_windows', '_state')

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


  def is_gap_window(self):

    if self._windows[WindowType.WGA].state == WindowType.N_WIN or\
      self._windows[WindowType.NO_WGA].state == WindowType.N_WIN:
        return True

    return False

  def get_window(self, wtype):
    return self._windows[wtype]


  """
  def get_rd_observations(self, name):

     if name == WindowType.BOTH:
          return (self._windows[WindowType.WGA].get_rd_observations(),
                  self._windows[WindowType.NO_WGA].get_rd_observations())
     elif name == WindowType.WGA:
          return self._windows[WindowType.WGA].get_rd_observations()
     elif name == WindowType.NO_WGA:
          return self._windows[WindowType.NO_WGA].get_rd_observations()

     raise Error("Name {0} is invalid ".format(name))
  """


  def get_rd_statistic(self, statistics="all", name=WindowType.BOTH):
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
          return (self._windows[WindowType.WGA].get_rd_statistic(statistic=statistics),
                  self._windows[WindowType.NO_WGA].get_rd_statistic(statistic=statistics))
        elif name == WindowType.WGA:
          return self._windows[WindowType.WGA].get_rd_statistic(statistic=statistics)
        elif name == WindowType.NO_WGA:
          return self._windows[WindowType.NO_WGA].get_rd_statistic(statistic=statistics)

        raise Error("Windowtype {0}"
                    " not in {1}".format(name, [WindowType.BOTH.name,
                                                WindowType.WGA.name,
                                                WindowType.NO_WGA.name]))











