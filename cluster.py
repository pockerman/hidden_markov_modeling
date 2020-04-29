"""
Cluster: A collection of windows
"""
from helpers import WindowState
from helpers import flat_windows_rd_from_indexes
from helpers import MixedWindowView
from preprocess_utils import compute_statistic



class Cluster(object):

  def __init__(self, id_, indexes):
    self._id = id_
    self._indexes = indexes
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

  def merge(self, cluster):
    self._indexes.extend(cluster.indexes)

  def get_data_from_windows(self, windows):
    return flat_windows_rd_from_indexes(indexes=self._indexes,
                                        windows=windows)

  def get_statistics(self, windows, statistic, window_type, **kwargs):

      wga_windows = [window.get_window(window_type)
                     for window in windows]

      window_data = self.get_data_from_windows(windows=wga_windows)
      return compute_statistic(data=window_data,statistics=statistic)

  def get_window_statistics(self, windows, statistic, **kwargs):
    statistics = []

    for idx in self.indexes:
      window = windows[idx]

      if isinstance(window, MixedWindowView):
        statistics.append(window.get_rd_stats(name=kwargs["window_type"],
                                              statistics=statistic))
      else:
        statistics.append(window.get_rd_stats(statistics=statistic))
    return statistics

  def get_bases(self, windows, windowtype="both"):

    bases = []
    for idx in self.indexes:
      window = windows[idx]
      bases.append(window.get_bases(windowtype=windowtype))
    return bases

