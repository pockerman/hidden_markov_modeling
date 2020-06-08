import argparse
import logging
import time
import json
import numpy as np
from multiprocessing import Process
from multiprocessing import Manager
from pomegranate import*

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import WindowType
from helpers import WindowState
from helpers import timefn
from helpers import INFO
from helpers import WARNING

from hmm_helpers import HMMCallback
from hmm_helpers import save_hmm

from region import Region

from cluster import Cluster
from cluster_utils import build_cluster_densities

from preprocess_utils import get_distributions_list_from_names as get_dist_list
from exceptions import Error

@timefn
def load_clusters(configuration):

  print("{0} Load clusters".format(INFO))
  clusters=[]

  for clst in configuration["clusters"]:
    cluster = Cluster.load(filename=configuration["clusters"][clst]["filename"])
    cluster.state = WindowState.from_string(configuration["clusters"][clst]["state"])
    clusters.append(cluster)

  return clusters

@timefn
def load_regions(configuration):

  print("{0} Load regions...".format(INFO))
  regions=[]

  for file in configuration["regions_files"]:
    region = Region.load(filename=file)

    if "check_windowing_sanity" in configuration and configuration["check_windowing_sanity"]:
        print("{0} Check window sanity for region {1}".format(INFO, region.ridx))
        region.check_windows_sanity()
        print("{0} Done...".format(INFO))
    region.get_mixed_windows()
    regions.append(region)

  return regions

@timefn
def init_hmm_2d(clusters, config):

  print("{0} Star HMM 2D initialization...".format(INFO))

  # get the HMM configuration
  hmm_config = config["HMM"]

  if WindowType.from_string(hmm_config["train_windowtype"]) !=\
        WindowType.BOTH:
          raise Error("For 2D Initialization WindowType shoulbe be 'BOTH'.")


  gap_state = None
  n_state_dist = None
  states = []


  if config["remove_windows_with_gaps"] == False:
    # we have a gap state then add it

    gap_windows_dist = config["gap_windows_dist"]
    n_state_dist = get_dist_list(dists_name=[gap_windows_dist["name"],
                                             gap_windows_dist["name"]],
                                       params={"uniform_params":gap_windows_dist["config"]["parameters"]})

    gap_state = \
            State(IndependentComponentsDistribution(n_state_dist), name="GAP_STATE")
    states.append(gap_state)


  # the rest of the states are MultivariateGaussians
  # we create as many states as clusters

  for cluster in clusters:

    if isinstance(cluster.state, str):
          name = cluster.state.lower()
    else:
          name = cluster.state.name.lower()

    kwargs={"exclude_gaps":True}
    mu1, mu2 = cluster.get_rd_statistics(statistic='mean',
                                         wtype=WindowType.BOTH,
                                         **kwargs)

    cov = cluster.get_rd_statistics(statistic='cov',
                                           wtype=WindowType.BOTH,
                                           **kwargs)
    state = State(MultivariateGaussianDistribution(means=np.array([mu1, mu2]),
                                                         covariance=cov), name=name)
    states.append(state)

  # print info about the states
  for state in states:
    print("{0} State: {1}".format(INFO, state.name))
    state_map = json.loads(str(state))
    print("{0} Distributions: {1}".format(INFO,
                                          state_map["distribution"]))

  # create the HMM
  hmm_model = HiddenMarkovModel(name=hmm_config["name"],
                                start=None, end=None)


  for state in states:
   prob = hmm_config["states"][state.name.lower()]["start_prob"]
   hmm_model.add_transition(hmm_model.start, state, prob)

  # add transitions for every state
  # to another this will create a dense HMM
  for i in states:
    for j in states:

      if i.name.lower()+"-"+j.name.lower() in hmm_config["transitions"]:
        prob =hmm_config["transitions"][i.name.lower()+"-"+j.name.lower()]
        hmm_model.add_transition(i, j, prob)
      else:
        print("{0} Transition from state"
              " {1} to state {2} is not specified".format(WARNING,
                                                          i.name.lower(),
                                                          j.name.lower()))

  # finally we need to bake
  hmm_model.bake(verbose=True)
  return hmm_model

@timefn
def init_hmm(clusters, configuration):


  print("{0} Star HMM initialization...".format(INFO))
  hmm_config = configuration["HMM"]
  n_state = None
  n_state_dist = None
  states = []

  # gap windows have not been removed
  if "remove_windows_with_gaps" in configuration and\
    configuration["remove_windows_with_gaps"] == False:

      # uniform distribution for gaps
      # so that E[X] = -999 and PDF = 1.0

      if WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.BOTH:

          n_state_dist = get_dist_list(dists_name=[configuration["gap_windows_dist"]["name"],
                                                   configuration["gap_windows_dist"]["name"]],
                                       params={"uniform_params":configuration["gap_windows_dist"]["config"]["parameters"]})

          n_state = \
            State(IndependentComponentsDistribution(n_state_dist), name="GAP_STATE")
      else:
        n_state_dist = get_dist_list(dists_name=[configuration["gap_windows_dist"]["name"],
                                                 configuration["gap_windows_dist"]["name"]],
                                     params={"uniform_params":configuration["gap_windows_dist"]["config"]["parameters"]})
        n_state = \
          State(n_state_dist[0], name="GAP_STATE")

  if n_state is not None:
    states.append(n_state)

  # create the HMM
  hmm_model = HiddenMarkovModel(name=hmm_config["name"],
                                start=None, end=None)

  for cluster in clusters:

    if isinstance(cluster.state, str):
          name = cluster.state.lower()
    else:
          name = cluster.state.name.lower()


    if WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.BOTH:
          states.append(State(IndependentComponentsDistribution([cluster.wga_density,
                                                                cluster.no_wga_density]),
                              name=name))
    elif WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.WGA:
          states.append(State(cluster.wga_density, name=name))
    elif WindowType.from_string(hmm_config["train_windowtype"]) ==\
        WindowType.NO_WGA:
          states.append(State(cluster.no_wga_density, name=name))
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
   prob = hmm_config["states"][state.name.lower()]["start_prob"]
   hmm_model.add_transition(hmm_model.start,
                            state, prob)

  # add transitions for every state
  # to another this will create a dense HMM
  for i in states:
    for j in states:

      if i.name.lower()+"-"+j.name.lower() in hmm_config["transitions"]:
        prob =hmm_config["transitions"][i.name.lower()+"-"+j.name.lower()]
        hmm_model.add_transition(i, j, prob)
      else:
        print("{0} Transition from state"
              " {1} to state {2} is not specified".format(WARNING,
                                                          i.name.lower(),
                                                          j.name.lower()))

  # finally we need to bake
  hmm_model.bake(verbose=True)
  return hmm_model

@timefn
def hmm_train(hmm_model, regions, configuration):

  print("{0} Star HMM training...".format(INFO))
  if configuration["HMM"]["train"] == True:
      print("{0} Creating training sequence...".format(INFO))

      hmm_conf = configuration["HMM"]
      if hmm_conf["train_sequence_source"] == "region":

        observations = []
        for region in regions:

          region_sequences = \
            region.get_region_as_rd_mean_sequences(size=hmm_conf["train_sequence_size"],
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
      print("{0} Fit HMM...".format(INFO))
      print("{0} Number of training sequences {1}".format(INFO,
                                                          len(observations)))

      print("{0} HMM transition matrix (before fit): ".format(INFO))
      print(hmm_model.dense_transition_matrix())
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
  hmm_model.bake(verbose=True)
  print("{0} Done...".format(INFO))

  if configuration["HMM"]["train"] == True:
    print("{0} HMM transition matrix (after fit): ".format(INFO))
  else:
    print("{0} HMM transition matrix: ".format(INFO))

  print(hmm_model.dense_transition_matrix())

  if configuration["HMM"]["save_model"]:
    print("{0} Saving HMM...".format(INFO))
    save_hmm(hmm_model=hmm_model,
                   configuration=configuration,
                   win_interval_length=0)
    print("{0} Done...".format(INFO))

  return hmm_model


def load_regions_worker(p, configuration, region_list,
                        msg_dict, errors_dict):

  try:
      region_list.extend(load_regions(configuration=configuration))
      msg_dict[p] = "Load regions worker finished"
  except Exception as e:
    msg = "An exception occured in load_region_worker. Exception string: " + str(e)
    errors_dict[p] = msg
    return

def load_clusters_worker(p, configuration, cluster_list,
                         msg_dict, errors_dict):

  try:
      cluster_list.extend(load_clusters(configuration=configuration))
      build_cluster_densities(clusters_lst=clusters, **configuration)
      msg_dict[p] = "Load clusters worker finished"
  except Exception as e:
    msg = "An exception occured in load_clusters_worker. Exception string: " + str(e)
    errors_dict[p] = msg
    return

@timefn
def main(configuration):

    print("{0} Set up logger".format(INFO))
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")
    print("{0} Done...".format(INFO))

    procs = []
    manager = Manager()
    regions_list = manager.list()
    errors_dict = manager.dict()
    errors_dict[0] = "No error"
    msg_dict = manager.dict()

    procs.append(Process(target=load_regions_worker,
                         args=(0, configuration, regions_list,
                               msg_dict, errors_dict)))
    procs[0].start()

    #regions_list = load_regions(configuration=configuration)

    # load the clusters whilst waiting
    # for the regions to load
    clusters = load_clusters(configuration=configuration)

    if configuration["HMM"]["use_multivariate"]:

      print("{0} Using Multivariate PDFs".format(INFO))

      # make sure the region has been loaded
      procs[0].join()

      if errors_dict[0] != "No error":
        raise Error(errors_dict[0])

      if len(regions_list) == 0:
        raise Error("Regions have not been loaded correctly")

      # assign the windows to the clusters
      for region in regions_list:
        region.get_mixed_windows()

      windows=[]

      for w in regions_list[0].get_mixed_windows():
        if w.is_gap_window():
          continue
        windows.append(w)

      for cluster in clusters:
        cluster.windows = windows

      hmm = init_hmm_2d(clusters=clusters,config=configuration)
    else:

      build_cluster_densities(clusters_lst=clusters, **configuration)
      hmm = init_hmm(clusters=clusters, configuration=configuration)

      # join here
      procs[0].join()

      if errors_dict[0] != "No error":
        raise Error(errors_dict[0])

      if len(regions_list) == 0:
        raise Error("Regions have not been loaded correctly")

    hmm = hmm_train(hmm_model=hmm,
                    regions=regions_list,
                    configuration=configuration)
    return hmm, regions_list


if __name__ == '__main__':
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

    print("{0} Read configuration file".format(INFO))
    args = parser.parse_args()
    configuration = read_configuration_file(args.config)
    print("{0} Done...".format(INFO))

    main(configuration=configuration)

    total_end = time.perf_counter()
    print("{0} Finished training. "
          "Execution time {1} secs".format(INFO, total_end - total_start))
