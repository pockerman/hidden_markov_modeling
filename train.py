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
from helpers import WindowType
from helpers import INFO
from region import Region
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

def make_window_regions(configuration):

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

        else:
            print("{0} No filtering windows"
                  " for Ns requested...".format(INFO))

        # compute the mixed windows for the region
        region.get_mixed_windows()

        print("{0} Number of mixed windows: {1}".format(INFO,
                                                        region.get_n_mixed_windows()))


        if "outlier_remove" in configuration and\
          configuration["outlier_remove"]:

            region.remove_outliers(configuration=configuration)
            print("{0} Number of windows after outlier removal: {1}".format(INFO,
                                                                        region.get_n_mixed_windows()))
        else:
          print("{0} No outlier removal performed".format(INFO))


        # save the region statistics
        region.save_mixed_windows_statistic(statistic="mean")
        regions_created.append(region)
        counter += 1

    return regions_created

def create_clusters(regions, configuration):

  kwargs = configuration["clusterer"]

  # assemble all the windows
  windows = []
  for region in regions:
    for window in region:
      windows.append(window)

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
                            indexes=clusters_indexes[i],
                            windows=windows))

  print("{0} Saving cluster indices".format(INFO))
  save_clusters(clusters=clusters, windows=windows, statistic="mean")
  print("{0} Done...".format(INFO))
  return clusters


def find_tuf_in_clusters(clusters, configuration):

  diff_rd = 0.0;
  tuf_cluster = None

  for cluster in clusters:
    mu1, mu2 = cluster.get_statistics(statistic="mean",
                                      window_type="both")

    if np.abs(mu1 - mu2) > diff_rd:
      diff_rd = np.abs(mu1 - mu2)
      tuf_cluster = cluster


  return diff_rd, tuf_cluster


def label_clusters(clusters, configuration):

  labeler = SignificanceTestLabeler(clusters=clusters, windows=windows)
  return labeler.label(test_config=configuration["labeler"])


def fit_clusters_distribution(clusters, configuration):

  kwargs = configuration["cluster_distribution"]
  print("{0} Fitting clusters densities...".format(INFO) )
  build_cluster_densities(clusters=clusters, **kwargs)
  print("{0} Done...".format(INFO))

  if configuration["save_cluster_densities"]:
    print("{0} Saving clusters densities...".format(INFO) )
    save_clusters_desnity(clusters=clusters, **kwargs)
    print("{0} Done...".format(INFO))


def init_hmm(clusters, configuration):

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

def hmm_train(clusters, regions, configuration):

  print("{0} HMM initialization...".format(INFO))

  hmm_model = init_hmm(clusters=clusters,
                       configuration=configuration)
  print("{0} Done...".format(INFO))
  print("{0} Creating training sequence...".format(INFO))

  if configuration["HMM"]["train_sequence_source"] == "region":

    observations = []
    for region in regions:
      observations.append(region.get_sequence(size=configuration["HMM"]["train_sequence_size"],
                                             window_type=WindowType.from_string(configuration["HMM"]["windowtype"])))

  elif configuration["HMM"]["train_sequence_source"] == "cluster":

    observations = []
    for cluster in clusters:
      observations.append(cluster.get_sequence(size=configuration["HMM"]["train_sequence_size"],
                                               window_type=WindowType.from_string(configuration["HMM"]["windowtype"])))

  else:
    raise Error("Training sequence type has not been specified")


  print("{0} Done...".format(INFO))

  print("{0} HMM transition matrix (before fit): ".format(INFO))
  print(hmm_model.dense_transition_matrix())

  print("{0} Fit HMM...".format(INFO))
  print("{0} Number of training sequences {1}".format(len(observations)))
  hmm_model, history = \
    hmm_model.fit(sequences=observations,
                  algorithm=configuration["HMM"]["train_solver"],
                  return_history=True,
                  verbose=configuration["HMM"]["verbose"],
                  lr_decay=configuration["HMM"]["lr_decay"],
                  callbacks=[HMMCallback(callback=print_logs_callback)],
                  inertia=configuration["HMM"]["inertia"])

  #finalize the model
  hmm_model.bake()
  print("{0} Done...".format(INFO))

  print("{0} HMM transition matrix (after fit): ".format(INFO))
  print(hmm_model.dense_transition_matrix())

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
    regions = make_window_regions(configuration=configuration)
    print("{0} Done...".format(INFO))

    print("{0} Start clustering....".format(INFO))
    clusters = create_clusters(regions=regions,
                               configuration=configuration)
    print("{0} Done...".format(INFO))

    print("{0} Compute max mean difference in clusters...".format(INFO))
    mean_diff, cluster = find_tuf_in_clusters(clusters=clusters,
                                              configuration=configuration)
    print("{0} Done...".format(INFO))
    print("{0} Max mean difference: {1} for cluster: {2} ".format(INFO,
                                                                  mean_diff,
                                                                  cluster.cidx))

    if configuration["label_clusters"]:
      print("{0} Labelling clusters...".format(INFO))
      clusters = label_clusters(clusters=clusters, windows=mixed_windows,
                                configuration=configuration)
      print("{0} Done...".format(INFO))

    print("{0} Fitting clusters distributions...".format(INFO))
    fit_clusters_distribution(clusters=clusters,
                              configuration=configuration)
    print("{0} Done...".format(INFO))
    print("{0} Starting HMM training...".format(INFO))
    hmm_train(clusters=clusters,
              regions=regions,
              configuration=configuration)

    print("{0} Done...".format(INFO))
    print("{0} Finished training".format(INFO))


if __name__ == '__main__':
    main()
