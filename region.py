from helpers import WindowType
from helpers import MixedWindowView
from exceptions import Error
from preprocess_utils import remove_outliers
from bam_helpers import extract_windows

class Region(object):
  """
  A region is simply a collection
  of windows
  """

  def __init__(self, idx, start, end, window_size):
    self._idx = idx
    self._start = start
    self._end = end
    self._w_size = window_size
    self._windows = {WindowType.WGA:[],
                     WindowType.NO_WGA:[] }

    self._mixed_windows = None


  def get_n_windows(self,type_):
    return len(self._windows[type_])

  def get_n_mixed_windows(self):
    if self._mixed_windows is None:
      return 0

    return len(self._mixed_windows)


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
    raise Error("Not implemented")

  def save_mixed_windows_statistic(self, statistic):
    raise Error("Not implemented")

  def remove_outliers(self, configuration):
    config = configuration["outlier_remove"]["config"]
    config["statistics"] = {"n_wga_w": no_wga_statistics,
                                  "wga_w":wga_statistics}

    self._mixed_windows = remove_outliers(windows=self._mixed_windows,
                                          removemethod=configuration["outlier_remove"]["name"],
                                          config=config)


