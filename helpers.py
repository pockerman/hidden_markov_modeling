import json
import numpy as np
import logging
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


Observation = namedtuple("Observation", ["position", "read_depth", "base"])


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
    def __init__(self, capacity):

        # maximum capacity of the window
        self._capacity = capacity

        # holds the observations i.e. the base
        # strings observed into the window
        self._observations = []

        # the total number of insertions/deletions
        self._net_indels = 0

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
        rd_data = [item[1] for item in self._observations]
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
        rd_data = [item[1] for item in self._observations]
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

    def insert_at(self, pos, data):
        """
        Insert the data at the specified position
        for the current window
        :param pos:
        :param data:
        :return:
        """
        self._observations.insert(pos, data)

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


