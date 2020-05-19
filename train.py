import argparse
import logging
import numpy as np
import time
import json
from pomegranate import*

from helpers import read_configuration_file
from helpers import set_up_logger

from helpers import Window
from helpers import WindowType
from helpers import WindowState
from helpers import INFO
from helpers import WARNING

from hmm_helpers import HMMCallback
from hmm_helpers import save_hmm

from region import Region
from analysis_helpers import save_clusters
from analysis_helpers import save_windows_statistic

from cluster import Cluster
from cluster_utils import build_cluster_densities, label_clusters
from preprocess_utils import build_clusterer
from preprocess_utils import get_distributions_list_from_names
from exceptions import Error

def load_clusters(configuration):
  clusters=[]

  for clst in configuration["clusters"]:
    cluster = Cluster.load(filename=clst["filename"])
    cluster.state = WindowState.from_string(clst["state"])
    clusters.append(cluster)

  return clusters


def init_hmm(clusters, configuration):


  hmm_config = configuration["HMM"]
  n_state = None
  n_state_dist = None
  states = []

  if "mark_N_windows" in configuration and\
    configuration["mark_N_windows"]:

      # uniform distribution for gaps
      # so that E[X] = -999 and PDF = 1.0

      if WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.BOTH:

          n_state_dist = get_distributions_list_from_names(dists_name=[configuration["n_windows_dist"]["name"],
                                                                      configuration["n_windows_dist"]["name"]],
                                                           params=configuration["n_windows_dist"]["config"]["parameters"])

          n_state = \
            State(IndependentComponentsDistribution(n_state_dist), name="GAP_STATE")
      else:
        n_state_dist = get_distributions_list_from_names(dists_name=[configuration["n_windows_dist"]["name"],
                                                                      configuration["n_windows_dist"]["name"]],
                                                           params=configuration["n_windows_dist"]["config"]["parameters"])
        n_state = \
          State(n_state_dist[0], name="GAP_STATE")

  if n_state is not None:
    states.append(n_state)

  # create the HMM
  hmm_model = HiddenMarkovModel(name=hmm_config["name"],
                                start=None, end=None)

  for cluster in clusters:
    name = cluster["name"]
    dists = cluster["distributions"]

    if WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.BOTH:
          states.append(State(IndependentComponentsDistribution(dists), name=name))
    elif WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.WGA:
          states.append(State(cluster.wga_density, name=use_name))
    elif WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.NO_WGA:
          states.append(State(cluster.no_wga_density, name=use_name))
    else:
      raise Error("Invalid train_windowtype. "
                  "{0} not in {1}".format(hmm_config["train_windowtype"],
                                          [WindowType.BOTH.name,
                                           WindowType.WGA.name,
                                           WindowType.NO_WGA.name]))

  for state in states:
    print("{0} State: {1}".format(INFO, state.name))
    state_map = json.loads(str(state))
    print("{0} Distributions: {1}".format(INFO,
                                          state_map["distribution"]))

  # add the states to the model
  hmm_model.add_states(states)

  for state in states:
   prob = hmm_config["states"][state.name]["start_prob"]
   hmm_model.add_transition(hmm_model.start,
                            state, prob)

  # add transitions for every state
  # to another this will create a dense HMM
  for i in states:
    for j in states:

      if i.name+"-"+j.name in hmm_config["transitions"]:
        prob =hmm_config["transitions"][i.name+"-"+j.name]
        hmm_model.add_transition(i, j, prob)
      else:
        print("{0} Transition from state"
              " {1} to state {2} is not specified".format(WARNING, i.name, j.name))

  # finally we need to bake
  hmm_model.bake(verbose=True)
  return hmm_model

def hmm_train(hmm_model, regions, configuration):


  if configuration["HMM"]["train"] == True:
      print("{0} Creating training sequence...".format(INFO))

      hmm_conf = configuration["HMM"]
      if hmm_conf["train_sequence_source"] == "region":

        observations = []
        for region in regions:

          region_sequences = \
            region.get_region_as_sequences(size=hmm_conf["train_sequence_size"],
                                           window_type=WindowType.from_string(hmm_conf["train_windowtype"]),
                                           n_seqs=hmm_conf["train_n_sequences_per_source"])

          for seq in region_sequences:
            observations.append(seq)

      elif hmm_conf["train_sequence_source"] == "cluster":

        observations = []
        for cluster in clusters:
          cluster_sequences = \
            cluster.get_sequence(size=hmm_conf["train_sequence_size"],
                                 window_type=WindowType.from_string(hmm_conf["train_windowtype"]),
                                 n_seqs=hmm_conf["train_n_sequences_per_source"])

          for seq in cluster_sequences:
            observations.append(seq)

      else:
        raise Error("Training sequence type has not been specified")

      print("{0} Done...".format(INFO))

      print("{0} HMM transition matrix (before fit): ".format(INFO))
      print(hmm_model.dense_transition_matrix())

      print("{0} Fit HMM...".format(INFO))
      print("{0} Number of training sequences {1}".format(INFO,
                                                          len(observations)))

      for i, seq in enumerate(observations):
        if -999.0 in seq:
          print("{0} Sequence {1} has -2.0".format(INFO, i))
          print(seq)
        elif (-999.0, -999.0) in seq:
          print("{0} Sequence {1} has -2.0".format(INFO, i))
          print(seq)

      print("{0} Training solver is: {1}".format(INFO,
                                                     configuration["HMM"]["train_solver"]))
      hmm_model, history = \
        hmm_model.fit(sequences=observations,
                      algorithm=configuration["HMM"]["train_solver"],
                      return_history=True,
                      verbose=configuration["HMM"]["verbose"],
                      lr_decay=configuration["HMM"]["lr_decay"],
                      callbacks=[HMMCallback(callback=None)],
                      inertia=configuration["HMM"]["inertia"])
  else:
      print("{0} No training performed.".format(INFO))

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

    print("{0} Start training...".format(INFO))
    total_start = time.perf_counter()
    description = "Check the README file for "
    "information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help="You must specify a json"
                        " formatted configuration file")


    print("{0} Read configuration file".format(INFO))
    args = parser.parse_args()
    configuration = read_configuration_file(args.config)
    print("{0} Done...".format(INFO))

    print("{0} Set up logger".format(INFO))
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")
    print("{0} Done...".format(INFO))

    print("{0} Load clusters".format(INFO))
    time_start = time.perf_counter()
    clusters = load_clusters(configuration=configuration)
    time_end = time.perf_counter()
    print("{0} Done. Execution time"
          " {1} secs".format(INFO, time_end - time_start))

    print("{0} Star HMM initialization...".format(INFO))
    time_start = time.perf_counter()

    hmm = init_hmm(clusters=clusters,configuration=configuration)
    time_end = time.perf_counter()
    print("{0} Done. Execution time"
          " {1} secs".format(INFO, time_end - time_start))

    print("{0} Star HMM training...".format(INFO))
    time_start = time.perf_counter()
    hmm_train(hmm_model=hmm,
              regions=regions,
              configuration=configuration)
    time_end = time.perf_counter()
    print("{0} Done. Execution time"
          " {1} secs".format(INFO, time_end - time_start))

    total_end = time.perf_counter()
    print("{0} Finished training. "
          "Execution time {1} secs".format(INFO, total_end - total_start))

if __name__ == '__main__':
    main()
