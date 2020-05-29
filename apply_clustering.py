import argparse
import logging
#import numpy as np
import time
import sys

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import WindowType
from helpers import INFO
from helpers import timefn

from region import Region
from analysis_helpers import save_clusters
#from analysis_helpers import save_windows_statistic
from analysis_helpers import save_clusters_gc_content

from cluster import Cluster
from cluster_utils import build_cluster_mean_and_std
from preprocess_utils import build_clusterer
from exceptions import Error

@timefn
def make_window_regions(configuration):

    print("{0} Creating window regions...".format(INFO))

    if configuration['processing']['type'] == 'multi':
      from parallel import par_make_window_regions
      return par_make_window_regions(configuration=configuration)

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
        sys.stdout.flush()
        print("{0} End index:   {1}".format(INFO, end_idx))
        sys.stdout.flush()
        region = Region(idx=counter,
                        start=start_idx,
                        end=end_idx,
                        window_size=windowsize)

        kwargs = {"sam_read_config":configuration["sam_read_config"]}

        if "debug" in configuration:
          kwargs["debug"] = configuration["debug"]

        print("{0} Creating WGA Windows...".format(INFO))
        sys.stdout.flush()
        region.make_wga_windows(chromosome=chromosome,
                                ref_filename=configuration["reference_file"]["filename"],
                                bam_filename=configuration["wga_file"]["filename"],
                                **kwargs)

        if region.get_n_windows(type_=WindowType.WGA) == 0:
            raise Error("WGA windows have not been created")
        else:
            print("{0} Number of WGA "
                  "windows: {1}".format(INFO,
                                        region.get_n_windows(type_=WindowType.WGA)))
            sys.stdout.flush()

        print("{0} Creating No WGA Windows...".format(INFO))
        sys.stdout.flush()
        region.make_no_wga_windows(chromosome=chromosome,
                                   ref_filename=configuration["reference_file"]["filename"],
                                   bam_filename=configuration["no_wga_file"]["filename"],
                                   **kwargs)

        if region.get_n_windows(type_=WindowType.NO_WGA) == 0:
            raise Error("Non-WGA windows have not  been created")
        else:
            print("{0} Number of Non WGA"
                  " windows: {1}".format(INFO,
                                         region.get_n_windows(type_=WindowType.NO_WGA)))
            sys.stdout.flush()

        regions_created.append(region)
        counter += 1

    return regions_created


def remove_gaps(region, configuration):
  # filter the windows for GAPs
    if "remove_windows_with_gaps" in configuration and\
          configuration["remove_windows_with_gaps"]:

            print("{0} Removing windows with gaps...".format(INFO))

            region.remove_windows_with_gaps()

            print("{0} Number of wga windows"
                  " after removing GAP windows: {1}".format(INFO,
                                                 region.get_n_windows(type_=WindowType.WGA)))
            sys.stdout.flush()
            print("{0} Number of non-wga windows"
                  " after removing GAP windows: {1}".format(INFO,
                                                 region.get_n_windows(type_=WindowType.NO_WGA)))
            print("{0} Done...".format(INFO))
            sys.stdout.flush()
    else:
            # mark the Gap windows
            print("{0} Marking Gap "
                  " windows with: {1}".format(INFO,
                                              configuration["mark_for_gap_windows"]))
            sys.stdout.flush()

            counter_ns = \
              region.mark_windows_with_gaps(n_mark=configuration["mark_for_gap_windows"])

            print("{0} Marked as Gap {1} Windows".format(INFO, counter_ns))
            sys.stdout.flush()


def remove_outliers(region, configuration):

  if "outlier_remove" in configuration and\
          configuration["outlier_remove"]:

            region.remove_outliers(configuration=configuration)
            print("{0} Number of windows "
                  "after outlier removal: {1}".format(INFO,
                                                      region.get_n_mixed_windows()))
            sys.stdout.flush()

            print("{0} Number of N windows "
                  "after outlier removal {1}".format(INFO,
                                                     region.count_gap_windows()))
            sys.stdout.flush()
  else:
          print("{0} No outlier "
                "removal performed".format(INFO))
          sys.stdout.flush()

@timefn
def clean_up_regions(regions, configuration):

  print("{0} Clean up regions".format(INFO))
  sys.stdout.flush()
  for region in regions:

    if "check_windowing_sanity" in configuration and \
          configuration["check_windowing_sanity"]:
            region.check_windows_sanity()

    # compute the mixed windows for the region
    region.get_mixed_windows()
    remove_gaps(region=region, configuration=configuration)

    print("{0} Number of mixed "
              "windows: {1}".format(INFO,
                                    region.get_n_mixed_windows()))
    sys.stdout.flush()

    print("{0} Number of GAP windows: {1}".format(INFO,
                                                    region.count_gap_windows()))
    sys.stdout.flush()

    remove_outliers(region=region, configuration=configuration)


def accumulate_windows(regions):
  # assemble all the windows
  windows = []
  for region in regions:
    for window in region:
      if not window.is_gap_window():
        windows.append(window)
  return windows


@timefn
def create_clusters(regions, configuration):

  print("{0} Start clustering....".format(INFO))
  sys.stdout.flush()
  kwargs = configuration["clusterer"]

  windows = accumulate_windows(regions=regions)

  # create the clusters
  clusterer, initial_index_medoids = \
    build_clusterer(data=windows,
                    nclusters=kwargs["config"]["n_clusters"],
                    method="kmedoids", **kwargs)

  print("{0} Initial medoids indexes: {1}".format(INFO,
                                                  initial_index_medoids))
  sys.stdout.flush()

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
  sys.stdout.flush()
  save_clusters(clusters=clusters, statistic="mean")

  if 'gc' in kwargs["config"]['features']:
    save_clusters_gc_content(clusters=clusters)

  print("{0} Done...".format(INFO))
  sys.stdout.flush()
  return clusters

@timefn
def make_clusters_mean_and_std(clusters, configuration):

  print("{0} Create clusters mean/std...".format(INFO))
  sys.stdout.flush()
  kwargs = {}
  print("{0} Make clusters mu/std...".format(INFO) )
  sys.stdout.flush()
  build_cluster_mean_and_std(clusters=clusters, **kwargs)
  print("{0} Done...".format(INFO))
  sys.stdout.flush()

def main(configuration):

    print("{0} Set up logger".format(INFO))
    sys.stdout.flush()
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")
    print("{0} Done...".format(INFO))
    sys.stdout.flush()

    regions = make_window_regions(configuration=configuration)

    clean_up_regions(regions=regions, configuration=configuration)

    print("{0} Saving regions...".format(INFO))
    sys.stdout.flush()
    time_start = time.perf_counter()

    for region in regions:
      region.save_mixed_windows_statistic(statistic="mean")
      region.save()
    time_end = time.perf_counter()
    print("{0} Done. Execution time {1} secs".format(INFO, time_end - time_start))
    sys.stdout.flush()

    clusters = create_clusters(regions=regions,
                               configuration=configuration)

    make_clusters_mean_and_std(clusters=clusters,
                               configuration=configuration)

    print("{0} Save clusters...".format(INFO))
    sys.stdout.flush()
    time_start = time.perf_counter()

    if "save_cluster_dbi" in configuration and\
      configuration["save_cluster_dbi"]:
        for cluster in clusters:
          cluster.diameter

          for other in clusters:
            cluster.calculate_distance_from_other(other=other)

    for cluster in clusters:
      cluster.save()
    time_end = time.perf_counter()
    print("{0} Done. Execution time {1} secs".format(INFO, time_end - time_start))
    sys.stdout.flush()

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
    sys.stdout.flush()

    main(configuration=configuration)
    total_end = time.perf_counter()
    print("{0} Finished clustering. "
          "Execution time {1} secs".format(INFO, total_end - total_start))
    sys.stdout.flush()
