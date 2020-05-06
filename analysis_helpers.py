"""
Various helpers to be used in the analysis module
"""

from helpers import WindowType

def save_cluster(filename, cluster, windows, statistic, wtype):

  with open(filename, 'w') as file:
    file.write(str(cluster.get_window_statistics(windows=windows,
                                                 statistic=statistic,
                                                 window_type=wtype)))

def save_clusters(clusters, windows, statistic):

  for cluster in clusters:
    wga_file = "cluster_"+str(cluster.cidx) +"_wga_w_" + statistic + ".txt"
    save_cluster(filename=wga_file, cluster=cluster,
                 windows=windows, statistic=statistic, wtype='wga_w')
    no_wga_file = "cluster_"+str(cluster.cidx) +"_no_wga_w_" + statistic + ".txt"
    save_cluster(filename=no_wga_file, cluster=cluster,
                 windows=windows, statistic=statistic, wtype='n_wga_w')

def save_windows_statistic(windows, statistic, region_id=None):

  window_stats = [window.get_rd_stats(statistics=statistic, name=WindowType.NO_WGA)
                        for window in windows]

  if region_id is not None:
    filename = "no_wga_windows_" + statistic + "_" + str(region_id) + ".txt"
  else:
    filename = "no_wga_windows_" + statistic + ".txt"

  with open(filename, 'w') as file:
    file.write(str(window_stats))

  window_stats = [window.get_rd_stats(statistics=statistic, name=WindowType.WGA)
                        for window in windows]

  if region_id is not None:
    filename = "wga_windows_" + statistic + "_" + str(region_id) + ".txt"
  else:
    filename = "wga_windows_" + statistic + ".txt"

  with open(filename, 'w') as file:
    file.write(str(window_stats))