"""
Cluster: A collection of windows
"""
from helpers import WindowState
from helpers import flat_windows_rd_from_indexes
from helpers import MixedWindowView
from helpers import WindowType
from helpers import WARNING
from preprocess_utils import compute_statistic



class Cluster(object):

  def __init__(self, id_, indexes, windows):
    self._id = id_
    self._indexes = indexes
    self._windows = windows
    self._state = WindowState.INVALID
    self._wga_density = None
    self._no_wga_density = None


  @property
  def state(self):
    return self._state

  @state.setter
  def state(self, value):
    self._state = value

  @property
  def cidx(self):
    return self._id

  @property
  def indexes(self):
    return self._indexes

  @property
  def wga_density(self):
    return self._wga_density

  @wga_density.setter
  def wga_density(self, value):
    self._wga_density = value

  @property
  def no_wga_density(self):
    return self._no_wga_density

  @no_wga_density.setter
  def no_wga_density(self, value):
    self._no_wga_density = value

  @property
  def windows(self):
    return self._windows

  def merge(self, cluster):
    self._indexes.extend(cluster.indexes)

  def get_data_from_windows(self):
    return flat_windows_rd_from_indexes(indexes=self._indexes,
                                        windows=self._windows)

  def get_sequence(self, size, window_type):

    sequence =[]

    if size < len(self._indexes):
        counter = 0
        for idx in self._indexes:
          window = self._windows[idx]
          sequence.append(window.get_rd_stats(statistics="mean"),
                          name=window_type)
          counter +=1

          if counter == size:
            break
    else:

      print("{0} Cluster size is less than {1}".format(WARNING, size))
      for idx in self._indexes:
          window = self._windows[idx]
          sequence.append(window.get_rd_stats(statistics="mean"),
                          name=window_type)

    return sequence

  def get_region_as_sequences(self, size, window_type, n_seqs):

    sequences = []
    sequence_local=[]
    for idx in self._indexes:
      window = self._windows[idx]
      sequence_local.append(window.get_rd_stats(statistics="mean", name=window_type))

      if len(sequence_local) == size:
        sequences.append(sequence_local)
        sequence_local=[]

      if n_seqs is not None and len(sequences) == n_seqs:
        break

    return sequences

  def get_statistics(self, statistic, window_type, **kwargs):


      if window_type == WindowType.BOTH:

        for index in self._indexes:
          window = self._windows[index]

          statistic1, statistic2 = \
            window.get_rd_stats(statistics=statistic)
          return statistic1, statistic2
      else:

        wga_windows = [window.get_window(window_type)
                       for window in self._windows]

        window_data = self.get_data_from_windows(windows=wga_windows)
        return compute_statistic(data=window_data,statistics=statistic)

  def get_window_statistics(self, statistic, **kwargs):
    statistics = []

    for idx in self.indexes:
      window = self._windows[idx]

      if isinstance(window, MixedWindowView):
        statistics.append(window.get_rd_stats(name=kwargs["window_type"],
                                              statistics=statistic))
      else:
        statistics.append(window.get_rd_stats(statistics=statistic))
    return statistics
