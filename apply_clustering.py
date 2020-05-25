import argparse
import logging
import numpy as np
import time
import sys

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import WindowType
from helpers import INFO
from helpers import timefn
from helpers import print_msg

from region import Region
from analysis_helpers import save_clusters
from analysis_helpers import save_windows_statistic

from cluster import Cluster
from cluster_utils import build_cluster_mean_and_std
from preprocess_utils import build_clusterer
from exceptions import Error

MASTER_PID = 0

def make_mpi_window_regions(configuration, pid, comm):

  #from mpi4py import MPI

  REGION_IDX_TAG = 11

  n_procs = comm.Get_size()

  if pid == 0:
    print_msg(msg="{0} Using {1} MPI processes: ".format(INFO, n_procs),
              pid=pid, master_pid=MASTER_PID)


    # check the regions. If more than one then this is
    # our partition. Otherwise if only one region then
    # we partition this region in P blocks
    regions = configuration["regions"]

    print_msg(msg="{0} Regions used {1}".format(INFO, regions),
              pid=pid, master_pid=MASTER_PID)

    regions_list = [ (start, end) for start, end
                    in zip(regions["start"], regions["end"])]

    if len(regions_list) == 1:
      # only one region. We partition this
      load = (regions_list[0][1] - regions_list[0][0])//n_procs

      print_msg(msg="{0} Load is: {1}".format(INFO, load),
              pid=pid, master_pid=MASTER_PID)

      # this one is working on
      this_start = regions_list[0][0]
      this_end = regions_list[0][0] + load

      print_msg(msg="{0} this_start is: {1}".format(INFO, this_start),
              pid=pid, master_pid=MASTER_PID)
      print_msg(msg="{0} this_end is: {1}".format(INFO, this_end),
              pid=pid, master_pid=MASTER_PID)

      # send to the other procs
      for p in range(0, n_procs-1):

        data = np.array([this_end + p*load, this_end + p*load + load], dtype=np.uint64)
        comm.send(data, dest=p+1, tag=REGION_IDX_TAG)

      print("{0} PID {1} got: {2}".format(INFO, p+1, data))
      print("{0} PID {1} Send to {2}".format(INFO, pid, p+1))

      # do the windowing for proc
      configuration["regions"]["start"]=[this_start]
      configuration["regions"]["end"]=[this_end]
      return make_window_regions(configuration=configuration,pid=pid)
  else:

      data = comm.recv(source=0, tag=REGION_IDX_TAG)
      #print("{0} PID {1} Received from {2}".format(INFO, pid, 0))
      configuration["regions"]["start"]=[data[0]]
      configuration["regions"]["end"]=[data[1]]
      return make_window_regions(configuration=configuration, pid=pid)


def make_window_regions(configuration, pid):

    print_msg(msg="{0} Creating window regions...".format(INFO),
              pid=pid, master_pid=MASTER_PID)

    windowsize = configuration["window_size"]
    chromosome = configuration["chromosome"]

    print_msg(msg="{0} PID {1} Window size: {2}".format(INFO, pid, windowsize),
              pid=pid, master_pid=MASTER_PID)

    print_msg(msg="{0} PID {1} Chromosome:  {2}".format(INFO, pid, chromosome),
              pid=pid, master_pid=MASTER_PID)

    regions = configuration["regions"]

    print_msg(msg="{0} PID {1} Regions used {2}".format(INFO, pid, regions),
              pid=pid, master_pid=MASTER_PID)

    regions_list = [ (start, end) for start, end
                    in zip(regions["start"], regions["end"])]

    regions_created = []

    counter=0
    for r in regions_list:

        start_idx = r[0]
        end_idx = r[1]

        print_msg(msg="{0} PID {1} Start index: {2}".format(INFO, pid, start_idx),
              pid=pid, master_pid=MASTER_PID)
        print_msg(msg="{0} PID {1} End index:   {2}".format(INFO, pid, end_idx),
              pid=pid, master_pid=MASTER_PID)

        region = Region(idx=counter,
                        pid=pid,
                        start=start_idx,
                        end=end_idx,
                        window_size=windowsize)

        kwargs = {"sam_read_config":configuration["sam_read_config"]}

        if "debug" in configuration:
          kwargs["debug"] = configuration["debug"]

        print_msg(msg="{0} PID {1} Creating WGA Windows...".format(INFO, pid),
              pid=pid, master_pid=MASTER_PID)

        region.make_wga_windows(chromosome=chromosome,
                                ref_filename=configuration["reference_file"]["filename"],
                                bam_filename=configuration["wga_file"]["filename"],
                                pid=pid,
                                **kwargs)

        if region.get_n_windows(type_=WindowType.WGA) == 0:
            raise Error("WGA windows have not been created")
        else:
            nwins = region.get_n_windows(type_=WindowType.WGA)
            print_msg(msg="{0} PID {1} Number of"
                      " WGA windows: {2}".format(INFO, pid, nwins),
                      pid=pid, master_pid=MASTER_PID)


        print_msg(msg="{0} PID {1} Creating No WGA Windows...".format(INFO, pid),
                  pid=pid, master_pid=MASTER_PID)

        region.make_no_wga_windows(chromosome=chromosome,
                                   ref_filename=configuration["reference_file"]["filename"],
                                   bam_filename=configuration["no_wga_file"]["filename"],
                                   pid=pid,
                                   **kwargs)

        if region.get_n_windows(type_=WindowType.NO_WGA) == 0:
            raise Error("Non-WGA windows have not  been created")
        else:
            nwins = region.get_n_windows(type_=WindowType.NO_WGA)

            print_msg(msg="{0} PID {1} Number of non-wga"
                  " windows: {2}".format(INFO, pid, nwins),
                  pid=pid, master_pid=MASTER_PID)

        if "check_windowing_sanity" in configuration and \
          configuration["check_windowing_sanity"]:

            region.check_windows_sanity()

        # compute the mixed windows for the region
        region.get_mixed_windows()

         # filter the windows for N's
        if "remove_windows_with_N" in configuration and\
          configuration["remove_windows_with_N"]:

            print("{0} Filtering windows for Ns...".format(INFO))

            region.remove_windows_with_ns()

            print("{0} Number of wga windows"
                  " after filtering: {1}".format(INFO,
                                                 region.get_n_windows(type_=WindowType.WGA)))
            print("{0} Number of non-wga windows"
                  " after filtering: {1}".format(INFO,
                                                 region.get_n_windows(type_=WindowType.NO_WGA)))
            print("{0} Done...".format(INFO))
        elif "mark_N_windows" in configuration and\
          configuration["mark_N_windows"]:

            print("{0} Marking N "
                  " windows with: {1}".format(INFO,
                                              configuration["mark_for_N_windows"]))
            counter_ns = \
              region.mark_windows_with_ns(n_mark=configuration["mark_for_N_windows"])

            print("{0} Marked as N {1} Windows".format(INFO, counter_ns))

        else:
            print("{0} No filtering windows"
                  " for Ns requested...".format(INFO))

        print("{0} Number of mixed "
              "windows: {1}".format(INFO,
                                    region.get_n_mixed_windows()))

        print("{0} Number of N windows: {1}".format(INFO,
                                                    region.count_n_windows()))

        regions_created.append(region)
        counter += 1

    return regions_created


def get_totals_from_regions(regions):

  total_sum = [0.0, 0.0]
  total_sqrd_sum = [0.0, 0.0]
  total_n_wins = 0
  for region in regions:

      total_n_wins += region.get_n_mixed_windows()
      rsum = region.get_rd_sum_statistic(wtype=WindowType.BOTH)
      total_sum[0] += rsum[0]
      total_sum[1] += rsum[1]
      rsqrdsum = region.get_rd_sum_sqrd_statistic(wtype=WindowType.BOTH)

      total_sqrd_sum[0] += rsqrdsum[0]
      total_sqrd_sum[1] += rsqrdsum[1]
  return total_sum, total_sqrd_sum, total_n_wins


def remove_outliers_from_regions(regions, configuration, pid, comm):

  if configuration["execution_type"] == "mpi":

    # we need to collect the global statistics
    total_sum, total_sqrd_sum, total_n_wins = get_totals_from_regions(regions=regions)

    # broadcast these to every participating
    # process
    recv = comm.alltoall((total_sum, total_sqrd_sum, total_n_wins))

    tsum = [total_sum[0], total_sum[1]]
    tsqrdsum = [total_sqrd_sum[0], total_sqrd_sum[1]]
    for item in recv:
      tsum[0] += item[0][0]
      tsum[1] += item[0][1]

      tsqrdsum[0] += item[1][0]
      tsqrdsum[1] += item[1][1]
      total_n_wins += item[2]

    total_mean = (tsum[0]/total_n_wins,
                  tsum[1]/total_n_wins)

    total_var = ((tsqrdsum[0] - tsum[0]**2)/total_n_wins,
                 (tsqrdsum[1] - tsum[1]**2)/total_n_wins)

    global_stats={"wga_statistics":{
        "mean": total_mean[0],
        "var":total_var[0]
      },
      "no_wga_statistics":{
          "mean": total_mean[1],
         "var":total_var[1]
        }
      }

  else:

    total_sum, total_sqrd_sum, total_n_wins = get_totals_from_regions(regions=regions)

    total_mean = (total_sum[0]/total_n_wins,
                  total_sum[1]/total_n_wins)

    total_var = ((total_sqrd_sum[0] - total_sum[0]**2)/total_n_wins,
                 (total_sqrd_sum[1] - total_sum[1]**2)/total_n_wins)

    global_stats={"wga_statistics":{
        "mean": total_mean[0],
        "var":total_var[0]
      },
      "no_wga_statistics":{
          "mean": total_mean[1],
         "var":total_var[1]
        }
      }

  for region in regions:
      region.remove_outliers(configuration=configuration, global_stats=global_stats)

      print_msg(msg="{0} Number of windows "
                  "after outlier removal: {1}".format(INFO,
                                                      region.get_n_mixed_windows()),
                  pid=pid, master_pid=MASTER_PID)

      print_msg(msg="{0} Number of N windows "
                  "after outlier removal {1}".format(INFO,
                                                     region.count_n_windows()),
                  pid=pid, master_pid=MASTER_PID)

@timefn
def create_clusters(regions, configuration):

  print("{0} Start clustering....".format(INFO))
  kwargs = configuration["clusterer"]

  # assemble all the windows
  windows = []
  for region in regions:
    for window in region:
      if not window.is_n_window():
        windows.append(window)

  # create the clusters
  clusterer, initial_index_medoids = \
    build_clusterer(data=windows,
                    nclusters=kwargs["config"]["n_clusters"],
                    method="kmedoids", **kwargs)

  print("{0} Initial medoids indexes: {1}".format(INFO,
                                                  initial_index_medoids))

  # get the window indexes
  clusters_indexes = clusterer.get_clusters()
  centers = clusterer.get_medoids()
  clusters = []

  for i in range(len(clusters_indexes)):
    clusters.append(Cluster(id_ = i,
                            indexes=clusters_indexes[i],
                            windows=windows,
                            center_idx=centers[i],
                            dist_metric=kwargs["config"]["metric"]))

  print("{0} Saving clusters means".format(INFO))
  save_clusters(clusters=clusters, statistic="mean")
  print("{0} Done...".format(INFO))
  return clusters

@timefn
def make_clusters_mean_and_std(clusters, configuration):

  print("{0} Create clusters mean/std...".format(INFO))
  kwargs = {}
  print("{0} Make clusters mu/std...".format(INFO) )
  build_cluster_mean_and_std(clusters=clusters, **kwargs)
  print("{0} Done...".format(INFO))

def main(configuration):


    pid = 0
    comm = None

    if configuration["execution_type"] == "mpi":
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      pid = comm.Get_rank()

    print_msg(msg="{0} PID {1} Set up logger".format(INFO, pid),
                  pid=pid, master_pid=MASTER_PID)

    set_up_logger(configuration=configuration, **{"pid":pid})
    logging.info("Checking if logger is sane...")
    print_msg(msg="{0} PID {1} Done...".format(INFO, pid),
                  pid=pid, master_pid=MASTER_PID)

    if configuration["execution_type"] == "mpi":
      regions = make_mpi_window_regions(configuration=configuration, pid=pid, comm=comm)
    else:
      regions = make_window_regions(configuration=configuration, pid=pid)

    if "outlier_remove" in configuration and\
          configuration["outlier_remove"]:
          remove_outliers_from_regions(regions=regions,
                                       configuration=configuration,
                                       pid=pid, comm=comm)
    else:
       print_msg(mg="{0} No outlier "
                "removal performed".format(INFO),
                  pid=pid, master_pid=MASTER_PID)


    print_msg(msg="{0} PID {1} Saving regions...".format(INFO, pid),
                  pid=pid, master_pid=MASTER_PID)

    time_start = time.perf_counter()

    #for region in regions:
    #  region.save_mixed_windows_statistic(statistic="mean")
    #  region.save()
    time_end = time.perf_counter()

    print_msg(msg="{0} Done. Execution time {1} secs".format(INFO, time_end - time_start),
                  pid=pid, master_pid=MASTER_PID)


    """
    clusters = create_clusters(regions=regions,
                               configuration=configuration)

    make_clusters_mean_and_std(clusters=clusters,
                               configuration=configuration)

    print("{0} Save clusters...".format(INFO))
    time_start = time.perf_counter()

    if "save_cluster_dbi" in configuration and\
      configuration["save_cluster_dbi"]:
        for cluster in clusters:
          cluster.diameter

          for other in clusters:
            cluster.distance_from_other(other=other)

    for cluster in clusters:
      cluster.save()
    """

    time_end = time.perf_counter()
    print_msg(msg="{0} Done. Execution time {1} secs".format(INFO, time_end - time_start),
                  pid=pid, master_pid=MASTER_PID)

if __name__ == '__main__':


    print("{0} Start clustering...".format(INFO))
    sys.stdout.flush()
    total_start = time.perf_counter()
    description = "Check the README file for "
    "information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help="You must specify a json "
                        "formatted configuration file")


    print("{0} Read configuration file".format(INFO))
    sys.stdout.flush()
    args = parser.parse_args()
    configuration = read_configuration_file(args.config)
    print("{0} Done...".format(INFO))

    main(configuration=configuration)
    total_end = time.perf_counter()
    print("{0} Finished clustering. "
          "Execution time {1} secs".format(INFO, total_end - total_start))
