"""
Unit tests for preprocess_utils module
"""

"""
Unit tests for helpers module
"""
import unittest
import json

from helpers import Window
from helpers import Observation
from helpers import WindowState
from preprocess_utils import z_score_window_clusterer
from preprocess_utils import ZScoreWindowCluster
from preprocess_utils import calculate_windows_zscore
from preprocess_utils import windows_tails_p


class TestPreprocessUtils(unittest.TestCase):


  def create_windows(self, nwindows, window_capacity):
    windows = []

    # create 5 windows
    counter = 0
    for i in range(nwindows):
      window = Window(idx=i, capacity=window_capacity)

      for obs in range(window.capacity()):

        observation = Observation(position=counter,
                                  read_depth=obs,
                                  base=['A'])

        window.add(observation=observation)
        counter += 1

      windows.append(window)
    return windows

  def test_calculate_windows_zscore(self):
      """
    Test if the windows score is calculated correctly

    """
      windows = self.create_windows(nwindows=5, window_capacity=5)
      self.assertEqual(len(windows), 5)
      zscores = calculate_windows_zscore(windows=windows)
      self.assertEqual(len(zscores), len(windows))

  def test_windows_tails_ps_1(self):

    windows = self.create_windows(nwindows=5, window_capacity=5)
    self.assertEqual(len(windows), 5)
    zscores = calculate_windows_zscore(windows=windows)
    self.assertEqual(len(zscores), len(windows))

    upper_ps, lower_ps = windows_tails_p(zscores=zscores, interval_length=1)
    self.assertEqual(len(upper_ps), len(zscores))
    self.assertEqual(len(lower_ps), len(zscores))

  def test_windows_tails_ps_2(self):

    windows = self.create_windows(nwindows=5, window_capacity=5)

    self.assertEqual(len(windows), 5)
    zscores = calculate_windows_zscore(windows=windows)
    self.assertEqual(len(zscores), len(windows))

    upper_ps, lower_ps = windows_tails_p(zscores=zscores, interval_length=2)
    self.assertEqual(len(upper_ps), 3)
    self.assertEqual(len(lower_ps), 3)

  def test_windows_tails_ps_3(self):

    windows = self.create_windows(nwindows=5,
                                  window_capacity=5)

    self.assertEqual(len(windows), 5)
    zscores = calculate_windows_zscore(windows=windows)
    self.assertEqual(len(zscores), len(windows))

    upper_ps, lower_ps = windows_tails_p(zscores=zscores,
                                         interval_length=3)

    self.assertEqual(len(upper_ps), 2)
    self.assertEqual(len(upper_ps[0]), 3)
    self.assertEqual(len(upper_ps[1]), 2)
    self.assertEqual(len(lower_ps), 2)
    self.assertEqual(len(lower_ps[0]), 3)
    self.assertEqual(len(lower_ps[1]), 2)

  def test_windows_tails_ps_4(self):

    windows = self.create_windows(nwindows=5, window_capacity=5)

    zscores = calculate_windows_zscore(windows=windows)
    self.assertEqual(len(zscores), len(windows))

    upper_ps, lower_ps = windows_tails_p(zscores=zscores,
                                         interval_length=4)
    self.assertEqual(len(upper_ps), 2)
    self.assertEqual(len(upper_ps[0]), 4)
    self.assertEqual(len(upper_ps[1]), 1)
    self.assertEqual(len(lower_ps), 2)
    self.assertEqual(len(lower_ps[0]), 4)
    self.assertEqual(len(lower_ps[1]), 1)


  def test_windows_tails_ps_5(self):

    windows = self.create_windows(nwindows=5, window_capacity=5)

    zscores = calculate_windows_zscore(windows=windows)
    self.assertEqual(len(zscores), len(windows))

    upper_ps, lower_ps = windows_tails_p(zscores=zscores,
                                         interval_length=5)
    self.assertEqual(len(upper_ps), 1)
    self.assertEqual(len(upper_ps[0]), 5)
    self.assertEqual(len(lower_ps), 1)
    self.assertEqual(len(lower_ps[0]), 5)

  def test_z_score_window_clusterer(self):
    """
    Test the ZScore clusterer. The test creates 5 windows with
    capacity of 5 observations. Each window then is filled with the
    following observations.

    Observation(position=counter, read_depth=0, base=['A'])
    Observation(position=counter, read_depth=1, base=['A'])
    Observation(position=counter, read_depth=2, base=['A'])
    Observation(position=counter, read_depth=3, base=['A'])
    Observation(position=counter, read_depth=4, base=['A'])

    We base statistics for the test only on read_depth so the
    position and counter are immateral. All windows should be
    clustered as NORMAL

    """

    clusterer = ZScoreWindowCluster(cutoff=0.05)

    # cretate the windows these will be of
    windows = self.create_windows(nwindows=5, window_capacity=5)
    self.assertEqual(len(windows), 5)
    windows = z_score_window_clusterer(windows=windows,
                                       n_consecutive_windows=1,
                                       selector=clusterer)


    self.assertEqual(len(windows), 5)
    for window in windows:
      self.assertEqual(window.get_state(), WindowState.NORMAL)



if __name__ == '__main__':
    unittest.main()