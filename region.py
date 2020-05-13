from helpers import WindowType
from helpers import MixedWindowView
from helpers import WARNING
from exceptions import Error
from preprocess_utils import remove_outliers
from preprocess_utils import compute_statistic
from analysis_helpers import save_windows_statistic
from bam_helpers import extract_windows


class RegionIterator(object):

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

class Region(object):
  """
  A region is simply a collection
  of windows
  """

  def __init__(self, idx, start, end, window_size):

    if end <= start:
      raise Error("Invalid start/end points. "
                  "End cannot be less than or equal to start")

    self._idx = idx
    self._start = start
    self._end = end
    self._w_size = window_size
    self._windows = {WindowType.WGA:[],
                     WindowType.NO_WGA:[] }

    self._mixed_windows = None


  def size(self):
    return self._end - self._start

  def get_n_windows(self,type_):
    return len(self._windows[type_])

  def get_n_mixed_windows(self):
    if self._mixed_windows is None:
      return 0

    return len(self._mixed_windows)

  def count_n_windows(self):

    counter = 0
    for win in self._mixed_windows:
      if win.is_n_window():
        counter += 1

    return counter


  def make_wga_windows(self, chromosome,
                       ref_filename,
                       test_filename, **kwargs):

    args = {"start_idx": self._start,
            "end_idx": self._end,
            "windowsize": self._w_size}

    if "quality_theshold" in kwargs:
      args["quality_theshold"] = kwargs["quality_theshold"]

    windows = extract_windows(chromosome=chromosome,
                                  ref_filename=ref_filename,
                                  test_filename=test_filename,
                                      **args)

    self._windows[WindowType.WGA] = windows

  def make_no_wga_windows(self, chromosome,
                          ref_filename,
                          test_filename, **kwargs):

    args = {"start_idx": self._start,
            "end_idx": self._end,
            "windowsize": self._w_size}

    if "quality_theshold" in kwargs:
      args["quality_theshold"] = kwargs["quality_theshold"]

    windows = extract_windows(chromosome=chromosome,
                                  ref_filename=ref_filename,
                                  test_filename=test_filename,
                                      **args)

    self._windows[WindowType.NO_WGA] = windows



  def get_mixed_windows(self):

    if self._mixed_windows is not None:
      return self._mixed_windows

    self._mixed_windows = []
    for win1, win2 in zip(self._windows[WindowType.WGA],
                          self._windows[WindowType.NO_WGA]):
          self._mixed_windows.append(MixedWindowView(wga_w=win1,
                                                     n_wga_w=win2))

    return self._mixed_windows

  def remove_windows_with_ns(self):

     # filter the windows for N's
     wga_filter_windows = [window
                           for window in self._windows[WindowType.WGA]
                                  if not window.has_base("N")]

     no_wga_filter_windows = [window
                              for window in self._windows[WindowType.NO_WGA]
                                  if not window.has_base("N")]

     self._windows[WindowType.WGA] = wga_filter_windows
     self._windows[WindowType.NO_WGA] = no_wga_filter_windows

  def save_mixed_windows_statistic(self, statistic):

    if self._mixed_windows is None:
      raise Error("Mixed windows have not been computed")

    save_windows_statistic(windows=self._mixed_windows,
                           statistic="mean", region_id=self._idx)

  def remove_outliers(self, configuration):

    if self._mixed_windows is None:
      raise Error("Mixed windows have not been computed")

    # compute the statistis

    wga_rds = []
    no_wga_rds = []

    for window in self._mixed_windows:
          if not window.is_n_window():
            wga_rds.extend(window.get_rd_counts(name=WindowType.WGA))
            no_wga_rds.extend(window.get_rd_counts(name=WindowType.NO_WGA))

    wga_statistics = compute_statistic(data=wga_rds,
                                       statistics="all")
    no_wga_statistics = compute_statistic(data=no_wga_rds,
                                          statistics="all")

    config = configuration["outlier_remove"]["config"]
    config["statistics"] = {WindowType.NO_WGA: no_wga_statistics,
                            WindowType.WGA:wga_statistics}

    self._mixed_windows = \
      remove_outliers(windows=self._mixed_windows,
                      removemethod=configuration["outlier_remove"]["name"],
                      config=config)

  def mark_windows_with_ns(self, n_mark):

    if self._mixed_windows is None:
      raise Error("Mixed windows have not been computed")

    counter=0
    for window in self._mixed_windows:
      wga_w = window.get_window(wtype=WindowType.WGA)
      n_wga_w = window.get_window(wtype=WindowType.NO_WGA)
      if wga_w.has_base("N") or n_wga_w.has_base("N"):
        wga_w.set_window_rd_mark(mark=n_mark)
        wga_w.state = WindowType.N_WIN

        n_wga_w.set_window_rd_mark(mark=n_mark)
        n_wga_w.state = WindowType.N_WIN
        counter += 1

    return counter

  def get_sequence(self, size, window_type):


    sequence =[]

    if size < len(self._mixed_windows):
        counter = 0
        for window in self._mixed_windows:
          sequence.append(window.get_rd_stats(statistics="mean", name=window_type))
          counter +=1

          if counter == size:
            break
    else:

      print("{0} Region size is less than {1}".format(WARNING, size))
      for window in self._mixed_windows:
          sequence.append(window.get_rd_stats(statistics="mean",name=window_type))

    return sequence

  def get_region_as_sequences(self, size, window_type, n_seqs):

    sequences = []
    sequence_local=[]
    for window in self._mixed_windows:
      sequence_local.append(window.get_rd_stats(statistics="mean",
                                                name=window_type))

      if len(sequence_local) == size:
        sequences.append(sequence_local)
        sequence_local=[]

      if n_seqs is not None and len(sequences) == n_seqs:
        break

    return sequences

  def __len__(self):
        return self.get_n_mixed_windows()

  def __iter__(self):
        """
        Produce an iterator to iterate over the accumulated
        window observations
        :return:
        """
        return RegionIterator(data=self._mixed_windows)

  def __getitem__(self, item):
        """
        Returns the item-th observation
        :param item: int the index of the observation to access
        :return: Observation instance
        """
        return self._mixed_windows[item]

  def __setitem__(self, o, value):
        """
        Set the o-th observation to value
        :param o:
        :param value:
        :return:
        """
        self._mixed_windows[o] = value


