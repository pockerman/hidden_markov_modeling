import argparse
import logging
import numpy as np
import time

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import WindowType
from helpers import INFO
from helpers import timefn

from region import Region
from analysis_helpers import save_clusters
from analysis_helpers import save_windows_statistic

from cluster import Cluster
from cluster_utils import build_cluster_mean_and_std
from preprocess_utils import build_clusterer
from exceptions import Error

@timefn
def make_window_regions(configuration):

    print("{0} Creating window regions...".format(INFO))

    windowsize = configuration["window_size"]
    chromosome = configuration["chromosome"]

    print("{0} Window size: {1}".format(INFO, windowsize))
    print("{0} Chromosome:  {1}".format(INFO, chromosome))

    regions = configuration["regions"]
    print("{0} Regions used {1}".format(INFO, regions))

    regions_list = [ (start, end) for start, end
                    in zip(regions["start"], regions["end"])]

    regions_created = []

    counter=0
    for r in regions_list:

        start_idx = r[0]
        end_idx = r[1]

        print("{0} Start index: {1}".format(INFO, start_idx))
        print("{0} End index:   {1}".format(INFO, end_idx))

        region = Region(idx=counter,
                        start=start_idx,
                        end=end_idx,
                        window_size=windowsize)

        kwargs = {}

        if "quality_theshold" in configuration:
          kwargs["quality_theshold"] = configuration["quality_theshold"]

        if "debug" in configuration:
          kwargs["debug"] = configuration["debug"]

        print("{0} Creating WGA Windows...".format(INFO))
        region.make_wga_windows(chromosome=chromosome,
                                ref_filename=configuration["reference_file"]["filename"],
                                test_filename=configuration["test_file"]["filename"],
                                **kwargs)

        if region.get_n_windows(type_=WindowType.WGA) == 0:
            raise Error("WGA windows have not been created")
        else:
            print("{0} Number of WGA windows: {1}".format(INFO,
                                                          region.get_n_windows(type_=WindowType.WGA)))

        print("{0} Creating No WGA Windows...".format(INFO))
        region.make_no_wga_windows(chromosome=chromosome,
                                   ref_filename=configuration["reference_file"]["filename"],
                                   test_filename=configuration["no_wga_file"]["filename"],
                                   **kwargs)

        if region.get_n_windows(type_=WindowType.NO_WGA) == 0:
            raise Error("Non-WGA windows have not  been created")
        else:
            print("{0} Number of non-wga"
                  " windows: {1}".format(INFO,
                                         region.get_n_windows(type_=WindowType.NO_WGA)))


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


        if "outlier_remove" in configuration and\
          configuration["outlier_remove"]:

            region.remove_outliers(configuration=configuration)
            print("{0} Number of windows "
                  "after outlier removal: {1}".format(INFO,
                                                      region.get_n_mixed_windows()))

            print("{0} Number of N windows "
                  "after outlier removal {1}".format(INFO,
                                                     region.count_n_windows()))

        else:
          print("{0} No outlier "
                "removal performed".format(INFO))


        # save the region statistics
        region.save_mixed_windows_statistic(statistic="mean")
        regions_created.append(region)
        counter += 1

    return regions_created

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


    print("{0} Set up logger".format(INFO))
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")
    print("{0} Done...".format(INFO))

    regions = make_window_regions(configuration=configuration)

    print("{0} Saving regions...".format(INFO))
    time_start = time.perf_counter()

    for region in regions:
      region.save()
    time_end = time.perf_counter()
    print("{0} Done. Execution time {1} secs".format(INFO, time_end - time_start))

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
    time_end = time.perf_counter()
    print("{0} Done. Execution time {1} secs".format(INFO, time_end - time_start))


if __name__ == '__main__':

    print("{0} Start clustering...".format(INFO))
    total_start = time.perf_counter()
    description = "Check the README file for "
    "information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help="You must specify a json "
                        "formatted configuration file")


    print("{0} Read configuration file".format(INFO))
    args = parser.parse_args()
    configuration = read_configuration_file(args.config)
    print("{0} Done...".format(INFO))

    main(configuration=configuration)
    total_end = time.perf_counter()
    print("{0} Finished clustering. "
          "Execution time {1} secs".format(INFO, total_end - total_start))
