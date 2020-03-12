import json
from collections import namedtuple
from exceptions import FullWindowException

def read_configuration_file(config_file):
    """
    Read the json configuration file and
    return a map with the configuration
    entries
    """
    with open(config_file) as json_file:
        configuration = json.load(json_file)
        return configuration


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

        self._capacity = capacity

        # holds the observations i.e. the base
        # strings observed into the window
        self._observations = []

        # the characteristic observation of the window
        self._mean_observation = None

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


