import argparse
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import csv
from pomegranate import*

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import save_hmm
from helpers import flat_windows
from helpers import flat_windows_from_state
from helpers import HMMCallback
from helpers import print_logs_callback
from helpers import flat_windows_rd_from_indexes
from helpers import MixedWindowView
from helpers import INFO
from analysis_helpers import save_clusters
from analysis_helpers import save_windows_statistic

from bam_helpers import extract_windows
from cluster import Cluster
from cluster_utils import clusters_statistics
from cluster_utils import build_cluster_densities
from cluster_utils import save_clusters_desnity
from hypothesis_testing import SignificanceTestLabeler
from preprocess_utils import fit_distribution
from preprocess_utils import compute_statistic
from preprocess_utils import build_clusterer
from preprocess_utils import remove_outliers
from exceptions import Error

def make_windows(configuration):

    wga_start_idx = configuration["test_file"]["start_idx"]
    wga_end_idx = configuration["test_file"]["end_idx"]

    if wga_end_idx == "none":
      wga_end_idx = None
    else:
      wga_end_idx = int(wga_end_idx)

    windowsize = configuration["window_size"]
    chromosome = configuration["chromosome"]

    print("{0} Start index: {1}".format(INFO, wga_start_idx))
    print("{0} End index:   {1}".format(INFO,wga_end_idx))
    print("{0} Window size: {1}".format(INFO, windowsize))
    print("{0} Chromosome:  {1}".format(INFO, chromosome))

    args = {"start_idx": int(wga_start_idx),
            "end_idx": wga_end_idx,
            "windowsize": int(windowsize)}

    if "quality_theshold" in configuration:
      args["quality_theshold"] = configuration["quality_theshold"]

    try:

        print("{0} Creating WGA Windows...".format(INFO))
        # extract the windows for the WGA treated file
        wga_windows = extract_windows(chromosome=chromosome,
                                      ref_filename=configuration["reference_file"]["filename"],
                                      test_filename=configuration["test_file"]["filename"],
                                      **args)

        if len(wga_windows) == 0:
            raise Error("WGA windows have not been created")
        else:
            print("{0} Number of WGA windows: {1}".format(INFO, len(wga_windows)))


        print("{0} Creating No WGA Windows...".format(INFO))
        non_wga_start_idx = configuration["no_wga_file"]["start_idx"]
        non_wga_end_idx = configuration["no_wga_file"]["end_idx"]

        args = {"start_idx": int(non_wga_start_idx),
                "end_idx": (non_wga_end_idx),
                "windowsize": int(windowsize)}

        if "quality_theshold" in configuration:
          args["quality_theshold"] = configuration["quality_theshold"]

        # exrtact the non-wga windows
        non_wga_windows = \
          extract_windows(chromosome=chromosome,
                          ref_filename=configuration["reference_file"]["filename"],
                          test_filename=configuration["no_wga_file"]["filename"],
                          **args)

        if len(non_wga_windows) == 0:
            raise Error("Non-WGA windows have not  been created")
        else:
            print("{0} Number of non-wga"
                  " windows: {1}".format(INFO,
                                         len(non_wga_windows)))


        # filter the windows for N's
        if "remove_windows_with_N" in configuration and\
          configuration["remove_windows_with_N"]:

            print("{0} Filtering windows for Ns...".format(INFO))
            wga_filter_windows = [window for window in wga_windows
                                  if not window.has_base("N")]

            no_wga_filter_windows = [window for window in non_wga_windows
                                  if not window.has_base("N")]

            wga_windows  = wga_filter_windows
            non_wga_windows = no_wga_filter_windows

            print("{0} Number of wga windows"
                  " after filtering: {1}".format(INFO,
                                                 len(wga_windows)))
            print("{0} Number of non-wga windows"
                  " after filtering: {1}".format(INFO,
                                                 len(non_wga_windows)))
            print("{0} Done...".format(INFO))

        else:
            print("{0} No filtering windows"
                  " for Ns requested...".format(INFO))


        # zip mixed windows the smallest length
        # prevails
        mixed_windows = []

        for win1, win2 in zip(wga_windows, non_wga_windows):
          mixed_windows.append(MixedWindowView(wga_w=win1,
                                               n_wga_w=win2))

        print("{0} Number of mixed windows: {1}".format(INFO,len(mixed_windows)))

        # compute the global statistics of the windows
        wga_rds = []
        no_wga_rds = []

        for window in mixed_windows:
          wga_rds.extend(window.get_rd_counts(name="wga_w"))
          no_wga_rds.extend(window.get_rd_counts(name="n_wga_w"))

        wga_statistics = compute_statistic(data=wga_rds, statistics="all")
        no_wga_statistics = compute_statistic(data=no_wga_rds, statistics="all")

        print("{0} WGA stats: {1}".format(INFO, wga_statistics))
        print("{0} No WGA stats: {1}".format(INFO, no_wga_statistics))

        save_windows_statistic(windows=mixed_windows, statistic="mean")


        # do the outlier removal

        if "outlier_remove" in configuration and\
          configuration["outlier_remove"]:

          config = configuration["outlier_remove"]["config"]
          config["statistics"] = {"n_wga_w": no_wga_statistics,
                                  "wga_w":wga_statistics}

          mixed_windows = remove_outliers(windows=mixed_windows,
                          removemethod=configuration["outlier_remove"]["name"],
                          config=config)

          print("{0} Number of windows after outlier removal: {1}".format(INFO,
                                                                          len(mixed_windows)))
        else:
          print("{0} No outlier removal performed".format(INFO))


        return mixed_windows

    except KeyError as e:
        logging.error("Key: {0} does not exit".format(str(e)))
        raise
    except Error as e:
        logging.error(str(e))
        raise
    except Exception as e:
        logging.error("Unknown exception occured: " + str(e))
        raise


def create_clusters(windows, configuration):

  kwargs = configuration["clusterer"]

  # create the clusters
  clusterer, initial_index_medoids = \
    build_clusterer(data=windows,
                    nclusters=len(configuration["states"]),
                    method="kmedoids", **kwargs)

  print("{0} Initial medoids indexes: {1}".format(INFO,
                                                  initial_index_medoids))

  # get the window indexes
  clusters_indexes = clusterer.get_clusters()
  clusters = []

  for i in range(len(clusters_indexes)):
    clusters.append(Cluster(id_ = i,
                            indexes=clusters_indexes[i]))

  print("{0} Saving cluster indices".format(INFO))
  save_clusters(clusters=clusters, windows=windows, statistic="mean")
  print("{0} Done...".format(INFO))
  return clusters


def label_clusters(clusters, windows, configuration):

  labeler = SignificanceTestLabeler(clusters=clusters, windows=windows)
  return labeler.label(test_config=configuration["labeler"])


def fit_clusters_distribution(clusters, windows, configuration):

  kwargs = configuration["cluster_distribution"]
  print("{0} Fitting clusters densities...".format(INFO) )
  build_cluster_densities(clusters=clusters,
                          windows=windows,
                          **kwargs)
  print("{0} Done...".format(INFO))

  if configuration["save_cluster_densities"]:
    print("{0} Saving clusters densities...".format(INFO) )
    save_clusters_desnity(clusters=clusters, windows=windows,
                          **kwargs)
    print("{0} Done...".format(INFO))


def init_hmm(clusters, windows, configuration):

  # create the HMM
  hmm_model = HiddenMarkovModel(name=configuration["HMM"]["name"],
                                start=None, end=None)


  state_to_dist = {}
  states = []
  for cluster, name in zip(clusters, configuration["states"]):
    states.append(State(IndependentComponentsDistribution([cluster.wga_density,
                                                           cluster.no_wga_density]),
                        name=name))


  # add the states to the model
  hmm_model.add_states(states)

  # construct the transition matrix.
  # We create a dense HMM. The probability
  # starting at a state is given in
  # configuration["HMM"]["start_prob"][state.name]

  for state in states:
      hmm_model.add_transition(hmm_model.start,
                             state,
                             configuration["HMM"]["start_prob"][state.name])

  # add transitions for every state
  # to another this will create a dense HMM
  for i in states:
    for j in states:

      if i.name == j.name:
        # high probabiity for self-transitioning
        hmm_model.add_transition(i, j, 0.95)
      else:

        #low probability for change state transition
        hmm_model.add_transition(i, j, 0.05)

  # finally we need to bake
  hmm_model.bake(verbose=True)
  return hmm_model

def hmm_train(clusters, windows, configuration):

  print("{0} HMM initialization...".format(INFO))

  hmm_model = init_hmm(clusters=clusters,
                       windows=windows,
                       configuration=configuration)
  print("{0} Done...".format(INFO))
  print("{0} Get observations from clusters...".format(INFO))

  observations = []
  windowtype = "both"

  for window in windows:

    if windowtype == "both":
      observations.append(window.get_rd_stats(statistics="mean",
                                              name=windowtype))
    else:
      observations.append([window.get_rd_stats(statistics="mean",
                                               name=windowtype)])

  print("{0} Done...".format(INFO))
  print("{0} Fit HMM...".format(INFO))
  hmm_model, history = \
    hmm_model.fit(sequences=observations,
                  algorithm=configuration["HMM"]["train_solver"],
                  return_history=True,
                  verbose=True,
                  lr_decay=0.6,
                  callbacks=[HMMCallback(callback=print_logs_callback)],
                  inertia=0.01)
  print("{0} Done...".format(INFO))

  if configuration["HMM"]["save_model"]:
    print("{0} Saving HMM...".format(INFO))
    save_hmm(hmm_model=hmm_model,
                   configuration=configuration,
                   win_interval_length=0)
    print("{0} Done...".format(INFO))


def main():

    print("{0} Starting training".format(INFO))
    description = "Check the README file for information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help='You must specify a json formatted configuration file')


    print("{0} Read configuratio file".format(INFO))
    args = parser.parse_args()
    configuration = read_configuration_file(args.config)
    print("{0} Done...".format(INFO))

    print("{0} Set up logger".format(INFO))
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")
    print("{0} Done...".format(INFO))

    print("{0} Creating windows...".format(INFO))
    mixed_windows = make_windows(configuration=configuration)
    print("{0} Done...".format(INFO))

    print("{0} Start clustering....".format(INFO))
    clusters = create_clusters(windows=mixed_windows,
                               configuration=configuration)
    print("{0} Done...".format(INFO))

    if configuration["label_clusters"]:
      print("{0} Labelling clusters...".format(INFO))
      clusters = label_clusters(clusters=clusters, windows=mixed_windows,
                                configuration=configuration)
      print("{0} Done...".format(INFO))

    print("{0} Fitting clusters distributions...".format(INFO))
    fit_clusters_distribution(clusters=clusters,
                              windows=mixed_windows,
                              configuration=configuration)
    print("{0} Done...".format(INFO))
    print("{0} Starting HMM training...".format(INFO))
    hmm_train(clusters=clusters,
              windows=mixed_windows,
              configuration=configuration)

    print("{0} Done...".format(INFO))
    print("{0} Finished training".format(INFO))


if __name__ == '__main__':
    main()
