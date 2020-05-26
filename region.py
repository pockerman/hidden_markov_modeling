from helpers import WindowType
from helpers import MixedWindowView
from helpers import Window
from helpers import WARNING, INFO
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

  @staticmethod
  def load(filename):
     with open(filename, 'r') as f:
       idx = int(f.readline().split(":")[1].rstrip("\n"))
       start = int(f.readline().split(":")[1].rstrip("\n"))
       end = int(f.readline().split(":")[1].rstrip("\n"))
       w_size = int(f.readline().split(":")[1].rstrip("\n"))

       region = Region(idx=idx, start=start,
                       end=end,window_size=w_size)

       n_wag_wins = int(f.readline().split(":")[1].rstrip("\n"))
       windows = []
       for w in range(n_wag_wins):
         wid = int(f.readline().split(":")[1].rstrip("\n"))
         cap = int(f.readline().split(":")[1].rstrip("\n"))
         size = int(f.readline().split(":")[1].rstrip("\n"))

         window = Window(idx=wid, capacity=cap)

         for obs in range(size):
           pos = int(f.readline().split(":")[1].rstrip("\n"))
           rd = float(f.readline().split(":")[1].rstrip("\n"))
           base = list(f.readline().split(":")[1].rstrip("\n"))

           #obs = Observation(position=pos, read_depth=rd, base=base)
           #window.add(observation=obs)
         windows.append(window)

       region.set_windows(wtype=WindowType.WGA, windows=windows)
       n_no_wag_wins = int(f.readline().split(":")[1].rstrip("\n"))

       windows = []
       for w in range(n_no_wag_wins):
         wid = int(f.readline().split(":")[1].rstrip("\n"))
         cap = int(f.readline().split(":")[1].rstrip("\n"))
         size = int(f.readline().split(":")[1].rstrip("\n"))

         window = Window(idx=wid, capacity=cap)

         for obs in range(size):
           pos = int(f.readline().split(":")[1].rstrip("\n"))
           rd = float(f.readline().split(":")[1].rstrip("\n"))
           base = list(f.readline().split(":")[1].rstrip("\n"))

           #obs = Observation(position=pos, read_depth=rd, base=base)
           #window.add(observation=obs)
         windows.append(window)
       region.set_windows(wtype=WindowType.NO_WGA, windows=windows)
       return region

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

  @property
  def ridx(self):
    return self._idx

  @property
  def w_size(self):
    return self._w_size

  @property
  def start(self):
    return self._start

  @property
  def end(self):
    return self._end

  def size(self):
    return self._end - self._start

  def get_n_windows(self,type_):
    return len(self._windows[type_])

  def get_n_mixed_windows(self):
    if self._mixed_windows is None:
      return 0

    return len(self._mixed_windows)

  def save(self):

    with open("region_" + str(self.ridx) + ".txt", 'w') as f:
      f.write("ID:"+str(self.ridx) + "\n")
      f.write("Start:"+str(self.start) + "\n")
      f.write("End:"+str(self.end) + "\n")
      f.write("WinSize:"+str(self.w_size) + "\n")

      f.write("WGA_N_WINDOWS:"+str(self.get_n_windows(type_=WindowType.WGA)) + "\n")

      for window in self._windows[WindowType.WGA]:
        f.write("WID:"+str(window.idx) + "\n")
        f.write("Capacity:"+str(window.capacity) + "\n")
        f.write("Size:"+str(len(window)) + "\n")
        f.write("N props:" +str(len(window.sam_property_names())) + "\n")
        for name in window.sam_property_names():
          f.write(name +":" + str(window.sam_property(name)) + "\n")

        f.write("N statistics" +str(len(window.get_statistics_map())) + "\n")

        for name in window.get_statistics_map():
          f.write(name +":" + str(window.get_statistics_value(name)) + "\n")

      f.write("NO_WGA_N_WINDOWS:"+str(self.get_n_windows(type_=WindowType.NO_WGA)) + "\n")
      for window in self._windows[WindowType.NO_WGA]:
        f.write("WID:"+str(window.idx) + "\n")
        f.write("Capacity:"+str(window.capacity) + "\n")
        f.write("Size:"+str(len(window)) + "\n")
        f.write("N props:" +str(len(window.sam_property_names())) +"\n")

        for name in window.sam_property_names():
          f.write(name +":" + str(window.sam_property(name)) + "\n")

        f.write("N statistics" +str(len(window.get_statistics_map())) + "\n")

        for name in window.get_statistics_map():
          f.write(name +":" + str(window.get_statistics_value(name)) + "\n")

  def count_n_windows(self):

    counter = 0
    for win in self._mixed_windows:
      if win.is_n_window():
        counter += 1

    return counter

  def set_windows(self, wtype, windows):

    if wtype != WindowType.WGA and\
      wtype != WindowType.NO_WGA:
        raise Error("Invalid Window type {0}"
                    " not in ['WGA', 'NO_WGA']".format(wtype))

    self._windows[wtype] = windows

  def make_wga_windows(self, chromosome,
                       ref_filename,
                       bam_filename, **kwargs):

    args = {"start_idx": self._start,
            "end_idx": self._end,
            "windowsize": self._w_size}

    if "debug" in kwargs:
      args["debug"] = kwargs["debug"]

    args["sam_read_config"]=kwargs["sam_read_config"]

    windows = extract_windows(chromosome=chromosome,
                              ref_filename=ref_filename,
                              bam_filename=bam_filename,
                                      **args)

    print("{0} Region Start Window Coords: Start/End idx {1}".format(INFO, windows[0].start_end_pos))
    print("{0} Region End Window Coords: Start/End idx {1}".format(INFO, windows[-1].start_end_pos))
    self._windows[WindowType.WGA] = windows

  def make_no_wga_windows(self, chromosome,
                          ref_filename,
                          bam_filename, **kwargs):

    args = {"start_idx": self._start,
            "end_idx": self._end,
            "windowsize": self._w_size}

    if "debug" in kwargs:
      args["debug"] = kwargs["debug"]

    args["sam_read_config"]=kwargs["sam_read_config"]

    windows = extract_windows(chromosome=chromosome,
                              ref_filename=ref_filename,
                              bam_filename=bam_filename,
                              **args)

    print("{0} Start Window: Start/End idx {1}".format(INFO, windows[0].start_end_pos))
    print("{0} End Window: Start/End idx {1}".format(INFO, windows[-1].start_end_pos))
    self._windows[WindowType.NO_WGA] = windows

  def check_windows_sanity(self):

    # check if the rest of the windows
    # are aligned
    self.get_mixed_windows()

    if len(self._windows[WindowType.NO_WGA]) > len(self._windows[WindowType.WGA]) :
      print("{0} Windows size mismatch"
            " WGA {1} NON_WGA {2}".format(WARNING,
                                         len(self._windows[WindowType.WGA]),
                                          len(self._windows[WindowType.NO_WGA])))
    elif len(self._windows[WindowType.NO_WGA]) < len(self._windows[WindowType.WGA]) :
        print("{0} Windows size mismatch"
            " WGA {1} NON_WGA {2}".format(WARNING,
                                          len(self._windows[WindowType.WGA]),
                                          len(self._windows[WindowType.NO_WGA])))


    for window in self._mixed_windows:
      start_wga, end_wga = window.get_window(wtype=WindowType.WGA).get_start_end_pos()
      start_no_wga, end_no_wga = window.get_window(wtype=WindowType.NO_WGA).get_start_end_pos()
      if (start_wga, end_wga) != (start_no_wga, end_no_wga):
          raise Error("Invalid window matching "
                      "window WGA at {0}, {1} "
                      "matched with NO WGA window at {2}, {3}".format(start_wga,
                                                                      end_wga,
                                                                      start_no_wga,
                                                                      end_no_wga))

  def get_mixed_windows(self):

    if self._mixed_windows is not None:
      return self._mixed_windows

    self._mixed_windows = []
    for win1, win2 in zip(self._windows[WindowType.WGA],
                          self._windows[WindowType.NO_WGA]):
          self._mixed_windows.append(MixedWindowView(wga_w=win1,
                                                     n_wga_w=win2))
    return self._mixed_windows

  def remove_windows_with_gaps(self):

     # filter the windows for N's
     wga_filter_windows = [window
                           for window in self._windows[WindowType.WGA]
                                  if not window.has_gaps()]

     no_wga_filter_windows = [window
                              for window in self._windows[WindowType.NO_WGA]
                                  if not window.has_gaps()]

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
            wga_rds.extend(window.get_rd_observations(name=WindowType.WGA))
            no_wga_rds.extend(window.get_rd_observations(name=WindowType.NO_WGA))


    if len(wga_rds) == 0 or len(no_wga_rds) == 0:
      print("{0} Cannot remove outliers for region. "
            "Empty RD list detected".format(WARNING))
      return

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

  def mark_windows_with_gaps(self, n_mark):

    if self._mixed_windows is None:
      raise Error("Mixed windows have not been computed")

    counter=0
    for window in self._mixed_windows:
      wga_w = window.get_window(wtype=WindowType.WGA)
      n_wga_w = window.get_window(wtype=WindowType.NO_WGA)

      ## Add error if one has N and the other not
      if wga_w.has_gaps() == True and\
        n_wga_w.has_gaps() == False:
        raise Error("WGA Window {0} has N "
                    "but Non WGA Window {1} does not".format(wga_w.idx,
                                                             n_wga_w.idx))
      elif wga_w.has_gaps() == False and\
        n_wga_w.has_gaps() == True:
        raise Error("WGA Window {0} does not have N "
                    "but Non WGA Window {1} does".format(wga_w.idx,
                                                         n_wga_w.idx))


      if wga_w.has_gaps() or n_wga_w.has_gaps():
        wga_w.set_window_rd_mark(mark=n_mark)
        wga_w.state = WindowType.N_WIN

        n_wga_w.set_window_rd_mark(mark=n_mark)
        n_wga_w.state = WindowType.N_WIN
        counter += 1

    return counter

  def get_rd_mean_sequence(self, size, window_type):


    sequence =[]

    if size < len(self._mixed_windows):
        counter = 0
        for window in self._mixed_windows:
          sequence.append(window.get_rd_statistic(statistics="mean", name=window_type))
          counter +=1

          if counter == size:
            break
    else:

      print("{0} Region size is less than {1}".format(WARNING, size))
      for window in self._mixed_windows:
          sequence.append(window.get_rd_statistic(statistics="mean",name=window_type))

    return sequence

  def get_region_as_rd_mean_sequences(self, size, window_type, n_seqs):

    sequences = []
    sequence_local=[]
    for window in self._mixed_windows:
      sequence_local.append(window.get_rd_statistic(statistics="mean",
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
        return RegionIterator(data=self._mixed_windows)

  def __getitem__(self, item):
        return self._mixed_windows[item]

  def __setitem__(self, o, value):
        self._mixed_windows[o] = value


