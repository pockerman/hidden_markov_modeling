import argparse
import logging
import time
import sys

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import WindowType
from helpers import INFO, WARNING
from helpers import timefn

from region import Region

from cluster import Cluster
from cluster_utils import build_cluster_mean_and_std
from preprocess_utils import build_clusterer
from exceptions import Error

@timefn
def make_window_regions(configuration):

    print("{0} Creating window regions...".format(INFO))
    print("{0} Processing type is: {1}".format(INFO, configuration['processing']['type']))

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


def remove_or_mark_gaps(region, configuration):
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
    else:
      print("{0} Window sanity check is disabled".format(WARNING))

    # compute the mixed windows for the region
    region.get_mixed_windows()
    remove_or_mark_gaps(region=region, configuration=configuration)

    print("{0} Number of mixed "
              "windows: {1}".format(INFO,
                                    region.get_n_mixed_windows()))
    sys.stdout.flush()

    print("{0} Number of GAP windows: {1}".format(INFO,
                                                    region.count_gap_windows()))
    sys.stdout.flush()

    remove_outliers(region=region, configuration=configuration)

@timefn
def save_regions(regions, configuration):

    print("{0} Saving regions".format(INFO))
    sys.stdout.flush()

    for region in regions:
      region.save_mixed_windows_statistic(statistic="mean")
      region.save()


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
                            dist_metric=kwargs["config"]["metric"],
                            dist_features=kwargs["config"]["features"]))

  print("{0} Saving clusters means".format(INFO))
  sys.stdout.flush()

  if configuration['processing']['type'] == 'multi':
    from multiprocessing import Process

    procs = []
    n_procs = configuration['processing']['n_procs']

    if len(clusters) < n_procs:
      for p in range(len(clusters)):
        procs.append(Process(target=save_clusters_content_worker,
                             args=(p, clusters[p], kwargs)))
        procs[p].start()
    elif len(clusters) == n_procs:
      for p in range(len(clusters) - 1):
        procs.append(Process(target=save_clusters_content_worker,
                             args=(p, clusters[p], kwargs)))
        procs[p].start()

      save_clusters_content_worker(p=len(clusters) - 1,
                                   cluster=clusters[len(clusters) - 1],
                                   kwrags=kwargs)
    else:

      for p in range(n_procs - 1):
        procs.append(Process(target=save_clusters_content_worker,
                             args=(p, clusters[p], kwargs)))
        procs[p].start()

      for c in range(n_procs - 1, len(clusters)):
        save_clusters_content_worker(p=c,
                                     cluster=clusters[c],
                                     kwrags=kwargs)

    for p in range(len(procs)):
       procs[p].join()

  else:
    from analysis_helpers import save_clusters

    save_clusters(clusters=clusters, statistic="mean",
                  tip=kwargs["config"]["metric"].upper())

    if 'gc' in kwargs["config"]['features']:
      print("{0} Saving clusters GC".format(INFO))
      sys.stdout.flush()
      from analysis_helpers import save_clusters_gc_content
      save_clusters_gc_content(clusters=clusters,
                               tip=kwargs["config"]["metric"].upper())

  return clusters

@timefn
def save_clusters_content_worker(p, cluster, kwargs):

  from analysis_helpers import save_cluster
  statistic='mean'
  tip = kwargs["config"]["metric"].upper()
  wga_file = "cluster_"+str(cluster.cidx) +"_wga_w_" + statistic + "_" + tip + ".txt"

  save_cluster(filename=wga_file,
               cluster=cluster,
               statistic=statistic, wtype=WindowType.WGA)

  no_wga_file = "cluster_"+str(cluster.cidx) +"_no_wga_w_" + statistic + "_" + tip + ".txt"
  save_cluster(filename=no_wga_file,
               cluster=cluster,
               statistic=statistic, wtype=WindowType.NO_WGA)

  if 'gc' in kwargs["config"]['features']:
    print("{0} Saving clusters GC".format(INFO))
    sys.stdout.flush()
    statistic = 'gc'
    wga_file = "cluster_"+str(cluster.cidx) +"_wga_w_" + statistic + "_" + tip + ".txt"
    save_cluster(filename=wga_file, cluster=cluster,
                 statistic=statistic, wtype=WindowType.WGA)


@timefn
def make_clusters_mean_and_std(clusters, configuration):

  print("{0} Create clusters mean/std...".format(INFO))
  sys.stdout.flush()
  kwargs = {}
  build_cluster_mean_and_std(clusters=clusters, **kwargs)


def save_cluster_worker(p, clusters, configuration):

  if "save_cluster_dbi" in configuration and\
      configuration["save_cluster_dbi"]:

          clusters[p].diameter

          for other in clusters:
            clusters[p].calculate_distance_from_other(other=other)

  clusters[p].save()


@timefn
def save_clusters(clusters, configuration):
   print("{0} Save clusters...".format(INFO))
   sys.stdout.flush()

   if configuration['processing']['type'] == 'multi':

     from multiprocessing import Process

     n_procs = configuration['processing']['n_procs']
     procs = []

     if len(clusters) < n_procs:
      for p in range(len(clusters)):
        procs.append(Process(target=save_cluster_worker,
                             args=(p, clusters, configuration)))
        procs[p].start()
     elif len(clusters) == n_procs:
      for p in range(len(clusters) - 1):
        procs.append(Process(target=save_cluster_worker,
                             args=(p, clusters, configuration)))
        procs[p].start()

      save_cluster_worker(p=len(clusters) - 1,
                          clusters=clusters,
                          configuration=configuration)
     else:

      for p in range(n_procs - 1):
        procs.append(Process(target=save_cluster_worker,
                             args=(p, clusters[p], configuration)))
        procs[p].start()

      for c in range(n_procs - 1, len(clusters)):
        save_cluster_worker(p=c, clusters=clusters, configuration=configuration)

     for p in range(len(procs)):
       procs[p].join()

     """
     for p in range(n_procs - 1):
       procs.append(Process(target=save_cluster_worker,
                            args=(p, clusters, configuration)))
       procs[p].start()
     save_cluster_worker(p=n_procs-1,
                         clusters=clusters,
                         configuration=configuration)

     for p in range(n_procs - 1):
       procs[p].join()
    """
   else:

     if "save_cluster_dbi" in configuration and\
        configuration["save_cluster_dbi"]:
          for cluster in clusters:
            cluster.diameter

            for other in clusters:
              cluster.calculate_distance_from_other(other=other)

     for cluster in clusters:
        cluster.save()


@timefn
def main(configuration):

    print("{0} Set up logger".format(INFO))
    sys.stdout.flush()
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")
    print("{0} Done...".format(INFO))
    sys.stdout.flush()

    regions = make_window_regions(configuration=configuration)

    clean_up_regions(regions=regions, configuration=configuration)

    proc = None
    if configuration['processing'] == 'multi':
      from multiprocessing import Process
      proc = Process(target=save_regions,
                     args=(regions, configuration))
      proc.start()
    else:
      save_regions(regions, configuration=configuration)

    clusters = create_clusters(regions=regions,
                               configuration=configuration)

    make_clusters_mean_and_std(clusters=clusters,
                               configuration=configuration)

    save_clusters(clusters=clusters,
                  configuration=configuration)

    if proc is not None:
      proc.join()


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
          "Total execution time {1} secs".format(INFO, total_end - total_start))
    sys.stdout.flush()
