"""
Module cluster
"""
import sys
import array
from pomegranate import*

from helpers import WindowState
from helpers import flat_windows_rd_from_indexes
from helpers import MixedWindowView
from helpers import WindowType
from helpers import WARNING, INFO
from exceptions import Error

from preprocess_utils import compute_statistic
from preprocess_utils import get_distance_metric



class Cluster(object):

  @staticmethod
  def load(filename):

    print("{0} Loading cluster from file: {1}".format(INFO, filename))

    with open(filename, 'r') as f:
      idx = int(f.readline().split(":")[1].rstrip("\n"))
      indices = f.readline().split(":")[1].rstrip("\n")
      indices = indices[1:-1].split(',')
      center = int(f.readline().split(":")[1].rstrip("\n"))
      dist_metric = f.readline().split(":")[1].rstrip("\n")
      dist_features = f.readline().split(":")[1].rstrip("\n")
      dist_features = dist_features.strip()
      dist_features = dist_features[1:-1].split(",")

      cluster = Cluster(id_=idx, indexes=[int(idx) for idx in indices],
                        windows=None, center_idx=center,
                        dist_metric=dist_metric.strip(),
                        dist_features=[ feature.strip().strip("'") for feature in dist_features])

      diameter = f.readline().split(":")[1].rstrip("\n")
      diameter = diameter.strip()
      if diameter != '-inf':
        diameter = float(diameter)
        cluster.diameter = diameter
      mean = float(f.readline().split(":")[1].rstrip("\n"))
      std = float(f.readline().split(":")[1].rstrip("\n"))
      cluster.wga_mean = mean
      cluster.wga_std = std

      mean = float(f.readline().split(":")[1].rstrip("\n"))
      std = float(f.readline().split(":")[1].rstrip("\n"))
      cluster.no_wga_mean = mean
      cluster.no_wga_std = std
      other_clusts = int(f.readline().split(":")[1].rstrip("\n"))

      if other_clusts >0 :
        dist_from_others={}
        for lidx in range(other_clusts):
          line = f.readline().split(":")
          indeces = line[0].split(",")
          idx1 = int(indeces[0])
          idx2 = int(indeces[1])
          value = float(line[1].rstrip("\n"))

          if idx1 != idx:
            raise Error("id {0} should have been "
                        "equal to cluster id {1} ".format(idx1, idx))

          dist_from_others[(idx1,idx2)]=value
        cluster.distance_from_others = dist_from_others

      return cluster

  @staticmethod
  def dbi(clusters):
    r={}

    for ci in clusters:
      for cj in clusters:

        if ci.cidx == cj.cidx:
          r[(ci.cidx, cj.cidx)] = 0.0
        else:
          si = ci.diameter
          sj = cj.diameter
          dij =  ci.get_distance_from_other(other=cj)
          r[(ci.cidx, cj.cidx)] = (si+sj)/dij

    # for each cluster find the maximum
    maxs = []
    for ci in clusters:
      max_ci = 0.0
      for cj in clusters:

        if r[(ci.cidx, cj.cidx)] > max_ci:
          max_ci = r[(ci.cidx, cj.cidx)]
      maxs.append(max_ci)

    return sum(maxs)/len(clusters)

  def __init__(self, id_, indexes, windows,
               center_idx, dist_metric, dist_features):
    self._id = id_
    self._indexes = indexes
    self._windows = windows
    self._center_idx = center_idx
    self._dist_metric = dist_metric
    self._dist_features = dist_features
    self._state = WindowState.INVALID
    self._wga_density = None
    self._wga_mean = 0.0
    self._wga_std = 0.0
    self._no_wga_density = None
    self._no_wga_mean = 0.0
    self._no_wga_std=0.0
    self._diameter = None
    self._distance_from_others=None

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
  def center_idx(self):
    return self._center_idx

  @property
  def center(self):
    return self._windows[self._center_idx]

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

  @no_wga_std.setter
  def no_wga_std(self, value):
    self._no_wga_std = value

  @property
  def windows(self):
    return self._windows

  @windows.setter
  def windows(self, value):
    self._windows = value


  @property
  def diameter(self):
    if self._diameter is not None:
      return self._diameter
    return self.calculate_diameter()

  @diameter.setter
  def diameter(self, value):
   self._diameter = value

  @property
  def distance_from_others(self):
    return self._distance_from_others

  @distance_from_others.setter
  def distance_from_others(self, value):
    self._distance_from_others = value

  def calculate_diameter(self):

    if self._windows is None:
      raise Error("Cannot calculate cluster "
                  "diameter without windows. "
                  "Cluster id is {0}".format(self.cidx))
    distance = 0.0

    for i in self._indexes:
      for j in self._indexes:
        w1 = self._windows[i]
        w2 = self._windows[j]

        if w1.is_gap_window() or w2.is_gap_window():
          continue

        w1_wga_mean, w1_no_wga_mean = w1.get_rd_statistic(statistics="mean",
                                                      name=WindowType.BOTH)

        w2_wga_mean, w2_no_wga_mean = w2.get_rd_statistic(statistics="mean",
                                                      name=WindowType.BOTH)

        metric = get_distance_metric(self._dist_metric, degree=4)
        new_dist = metric([w1_wga_mean, w1_no_wga_mean], [w2_wga_mean, w2_no_wga_mean])

        if new_dist  > distance:
          distance = new_dist
    self._diameter = distance
    return distance

  def get_distance_from_other(self, other):

    if self._distance_from_others is None:
      return "-inf"

    if (self.cidx, other.cidx) not in self._distance_from_others:
      return False

    return self._distance_from_others[(self.cidx, other.cidx)]

  def set_distance_from_other(self, other, dist):
    if self._distance_from_others is None:
        self._distance_from_others = {(self.cidx,other.cidx):dist}
    else:
        self._distance_from_others[(self.cidx,other.cidx)] = dist


  def calculate_distance_from_other(self, other):

    if other.cidx == self.cidx:
      self.set_distance_from_other(other, 0.0)
      return 0.0

    if self._distance_from_others is not None and \
      (self.cidx, other.cidx) in self._distance_from_others:
        return self._distance_from_others[(self.cidx, other.cidx)]

    if other.get_distance_from_other(other=self) != '-inf' or\
      other.get_distance_from_other(other=self) != False:
        self.set_distance_from_other(other,
                                     other.get_distance_from_other(other=self))
        return other.get_distance_from_other(other=self)

    distance = sys.float_info.max
    for i in self._indexes:
      for j in other.indexes:
        this_w = self._windows[i]
        other_w = self._windows[j]

        if this_w.is_gap_window() or other_w.is_gap_window():
          continue

        w1_wga_mean, w1_no_wga_mean = this_w.get_rd_statistic(statistics="mean",
                                                              name=WindowType.BOTH)

        w2_wga_mean, w2_no_wga_mean = other_w.get_rd_statistic(statistics="mean",
                                                               name=WindowType.BOTH)

        metric = get_distance_metric(self._dist_metric, degree=4)
        new_dist = metric([w1_wga_mean, w1_no_wga_mean], [w2_wga_mean, w2_no_wga_mean])

        if new_dist  < distance:
          distance = new_dist

    self.set_distance_from_other(other=other, dist=distance)
    return distance

  def distance_from(self, item):

    import copy

    cluster_center = self.center
    features = copy.deepcopy(self._dist_features)

    has_gc = False
    if 'gc' in features:
      features.pop(features.index('gc'))
      has_gc = True


    item_values = item.get_features(features=features)
    center_values = cluster_center.get_features(features=features)

    if has_gc:
      item_values.append(item.get_feature(feature='gc')[0])
      center_values.append(cluster_center.get_feature(feature='gc')[0])

    from preprocess_utils import get_distance_metric
    metric = get_distance_metric(dist_metric=self._dist_metric)
    distance = metric(item_values, center_values)
    return distance

  def merge(self, cluster):
    self._indexes.extend(cluster.indexes)

  def save(self, tips):
    filename = "cluster_" + str(self.cidx)
    if tips is not None:
      for tip in tips:
        filename += "_" + tip

    filename += ".txt"

    with open(filename, 'w') as f:

      f.write("ID:"+str(self.cidx) +"\n")
      f.write("Indices:" + str(self.indexes) +"\n")
      f.write("Center:"  + str(self._center_idx) +"\n")
      f.write("Dist Metric: " + self._dist_metric + "\n")
      f.write("Dist Features: " + str(self._dist_features) + "\n")

      diam =  str(self._diameter) if self._diameter is not None else '-inf'
      f.write("Diameter:" + diam + "\n")
      f.write("WGA_MEAN:"+str(self.wga_mean) +"\n")
      f.write("WGA_STD:" + str(self.wga_std) +"\n")
      f.write("NO_WGA_MEAN:" + str(self.no_wga_mean) +"\n")
      f.write("NO_WGA_STD:" + str(self.no_wga_std) +"\n")

      if self._distance_from_others is None:
        f.write("N_DIST_OTHERS:" + str(0) +"\n")
      else:
        f.write("N_DIST_OTHERS:" + str(len(self._distance_from_others)) +"\n")

        for item in self._distance_from_others:
          f.write(str(item[0]) + "," + str(item[1])+":" + str(self._distance_from_others[item]) +"\n")

  def get_sequence(self, size, window_type):

    sequence =[]

    if size < len(self._indexes):
        counter = 0
        for idx in self._indexes:
          window = self._windows[idx]
          sequence.append(window.get_rd_statistic(statistics="mean"),
                          name=window_type)
          counter +=1

          if counter == size:
            break
    else:

      print("{0} Cluster size is less than {1}".format(WARNING, size))
      for idx in self._indexes:
          window = self._windows[idx]
          sequence.append(window.get_rd_statistic(statistics="mean"),
                          name=window_type)

    return sequence

  def get_region_as_sequences(self, size, window_type, n_seqs):

    sequences = []
    sequence_local=[]
    for idx in self._indexes:
      window = self._windows[idx]
      sequence_local.append(window.get_rd_statistic(statistics="mean", name=window_type))

      if len(sequence_local) == size:
        sequences.append(sequence_local)
        sequence_local=[]

      if n_seqs is not None and len(sequences) == n_seqs:
        break

    return sequences


  def get_rd_sequnce(self, statistic, wtype, exclude_gaps):

    if self._indexes is None:
      raise Error("No indexes have been given to the cluster")

    if self._windows is None:
      raise Error("No windows have been given to the cluster")

    sequence = array.array('d')

    for index in self._indexes:
          window = self._windows[index]

          if exclude_gaps and window.is_gap_window():
                continue

          stat = window.get_rd_statistic(statistics=statistic, name=wtype)
          sequence.append(stat)
    return sequence


  def get_rd_statistics(self, statistic, wtype, **kwargs):

      if statistic == 'cov':

        if wtype != WindowType.BOTH:
          raise Error("Calculating covariance requires both samples")

        wdata_wga = self.get_rd_sequnce(statistic=statistic,
                                        wtype=WindowType.WGA,
                                        exclude_gaps=kwargs["exclude_gaps"])
        wdata_no_wga = self.get_rd_sequnce(statistic=statistic,
                                           wtype=WindowType.NO_WGA,
                                           exclude_gaps=kwargs["exclude_gaps"])

        X = np.stack((wdata_wga, wdata_no_wga), axis=0)
        return np.cov(X)


      if wtype == WindowType.BOTH:

        wdata_wga = self.get_rd_sequnce(statistic=statistic,
                                        wtype=WindowType.WGA,
                                        exclude_gaps=kwargs["exclude_gaps"])
        wdata_no_wga = self.get_rd_sequnce(statistic=statistic,
                                           wtype=WindowType.NO_WGA,
                                           exclude_gaps=kwargs["exclude_gaps"])

        stat1 = compute_statistic(data=wdata_wga, statistics=statistic)
        stat2 = compute_statistic(data=wdata_no_wga, statistics=statistic)
        return stat1, stat2

      raise Error("Invalid operation in cluster")

      """
      else:


        windows = [window.get_window(wtype)
                       for window in self._windows]

        window_data = flat_windows_rd_from_indexes(indexes=self._indexes,
                                                   windows=windows)

      return compute_statistic(data=window_data, statistics=statistic)
      """

  def get_window_statistics(self, statistic, **kwargs):


    statistics = []

    for idx in self.indexes:
      window = self._windows[idx]

      if isinstance(window, MixedWindowView):
        if statistic == 'gc':
          statistics.append(window.get_feature(feature='gc',name=kwargs["window_type"]))
        else:
          statistics.append(window.get_rd_statistic(name=kwargs["window_type"],
                                              statistics=statistic))
      else:
        statistics.append(window.get_rd_statistic(statistics=statistic))
    return statistics
