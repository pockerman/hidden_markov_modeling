"""
Cluster: A collection of windows
"""
from helpers import WindowState
from helpers import flat_windows_rd_from_indexes
from helpers import MixedWindowView
from pomegranate import*
from helpers import WindowType
from helpers import WARNING
from preprocess_utils import compute_statistic



class Cluster(object):

  @staticmethod
  def load(filename):
    with open(filename, 'r') as f:
      idx = int(f.readline().split(":")[1])
      indices = list(f.readline().split(":")[1])

      cluster = Cluster(id_=idx, indexes=indices, windows=None)
      mean = float(f.readline().split(":")[1])
      std = float(f.readline().split(":")[1])
      cluster.wga_mean = mean
      cluster.wga_std = std

      mean = float(f.readline().split(":")[1])
      std = float(f.readline().split(":")[1])
      cluster.no_wga_mean = mean
      cluster.no_wga_std = std

      return cluster

  def __init__(self, id_, indexes, windows):
    self._id = id_
    self._indexes = indexes
    self._windows = windows
    self._state = WindowState.INVALID
    self._wga_density = None
    self._wga_mean = 0.0
    self._wga_std = 0.0
    self._no_wga_density = None
    self._no_wga_mean = 0.0
    self._no_wga_std=0.0


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
  def wga_mean(self):
    return self._wga_mean

  @wga_mean.setter
  def wga_mean(self, value):
    self._wga_mean = value

  @property
  def wga_std(self):
    return self._wga_std

  @wga_std.setter
  def wga_std(self, value):
    self._wga_std = value

  @property
  def no_wga_density(self):
    return self._no_wga_density

  @no_wga_density.setter
  def no_wga_density(self, value):
    self._no_wga_density = value

  @property
  def no_wga_mean(self):
    return self._no_wga_mean

  @no_wga_mean.setter
  def no_wga_mean(self, value):
    self._no_wga_mean = value

  @property
  def no_wga_std(self):
    return self._no_wga_std

  @wga_std.setter
  def no_wga_std(self, value):
    self._no_wga_std = value

  @property
  def windows(self):
    return self._windows

  def merge(self, cluster):
    self._indexes.extend(cluster.indexes)


  def save(self):
    with open("cluster_" + str(self.cidx) + ".txt", 'w') as f:

      f.write("ID:"+str(self.cidx)+"\n")
      f.write("Indices:"+str(self.indexes)+"\n")
      f.write("WGA_MEAN:"+str(self.wga_mean)+"\n")
      f.write("WGA_STD:"+str(self.wga_std)+"\n")
      f.write("NO_WGA_MEAN:"+str(self.no_wga_mean)+"\n")
      f.write("NO_WGA_STD:"+str(self.no_wga_std)+"\n")


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

        window_data = flat_windows_rd_from_indexes(indexes=self._indexes,
                                                   windows=wga_windows)

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
