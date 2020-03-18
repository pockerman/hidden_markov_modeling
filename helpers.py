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
        return self.capacity() - len(self._observations) != 0

    def get_range(self, start, end):
        return self._observations[start:end]

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
        #elif at_count == 0:
        #    return 1

        return gc_count / (gc_count + at_count)

    def set_net_indels(self, net_indels):
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

    def insert_at_positions(self, data):
        """
        Insert the data
        :param data:
        :return:
        """
        for item in data:
            self.insert_at(pos=item[3] -1, data=item[:3])


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


