"""
Use synthetic data to evaluate implementation
and assumptions. The following states are
generated

NORMAL
DELETE
INSERT

"""

from pomegranate import*
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from helpers import read_configuration_file
from helpers import Window
from helpers import WindowState
from helpers import Observation
from helpers import add_window_observation
from helpers import DUMMY_ID
from helpers import windows_to_json
from helpers import flat_windows
from helpers import flat_windows_from_state
from helpers import HMMCallback
from helpers import print_logs_callback
from helpers import windows_rd_statistics
from exceptions import Error
from bam_helpers import DUMMY_BASE
from preprocess_utils import cluster
from preprocess_utils import fit_distribution


def create_synthetic_data(configuration, create_windows=True):

  bases = ['A', 'T', 'C', 'G']
  flip = ["YES", "NO"]
  states = {"normal_rd": configuration["normal_rd"],
            "delete_rd": configuration["delete_rd"],
            "insert_rd": configuration["insert_rd"]}

  seq_size = configuration["sequence_size"]

  synthetic_data = []
  state_name = np.random.choice(list(states.keys()))

  counter = 0
  for i in range(seq_size):

    if np.random.choice(flip) == 'YES' and \
      counter > configuration["minimum_occurences"]:
      state_name = np.random.choice(list(states.keys()))
      synthetic_data.append(states[state_name])
      counter = 0
      counter += 1
    #elif counter  configuration["swap_state_freq"]:

     # state_name = np.random.choice(list(states.keys()))
     # synthetic_data.append(states[state_name])
     # counter += 1

    else:
      synthetic_data.append(states[state_name])
      counter += 1


  if len(synthetic_data) != seq_size:
    raise Error("Invalid size of synthetic \
                data {0} not equal to {1}".format(len(synthetic_data), seq_size))

  if not create_windows:
    return synthetic_data

  windows = []

  idstart = 0
  window = Window(idx=idstart, capacity=configuration["window_size"])
  previous_observation = None

  for idx, item in enumerate(synthetic_data):

    # create an observation
    observation = Observation(position=idx,
                                  read_depth=item,
                                  base=bases)
    window = add_window_observation(window=window, windows=windows,
                                        observation=observation,
                                        windowcapacity=configuration["window_size"])


  if len(window) != window.capacity():

      # fill in missing data if this is requested
      if configuration.get("fill_missing_window_data", False):
        while window.has_capacity():
          window.add(observation=Observation(position=DUMMY_ID,
                                  read_depth=configuration["fill_missing_window_data_factor"],
                                  base= [DUMMY_BASE]))

  windows.append(window)
  return windows


def init_hmm(hmm_model, windows, configuration):

  # collect the windows with the same state
  normal_state = []
  delete_state = []
  insert_state = []
  states_to_windows = {}

  for window in windows:
        if window.get_state() == WindowState.NORMAL:
          normal_state.extend(window.get_rd_observations())
        elif window.get_state() == WindowState.DELETE:
          delete_state.extend(window.get_rd_observations())
        elif window.get_state() == WindowState.INSERT:
          insert_state.extend(window.get_rd_observations())

  states_to_windows["NORMAL"] = normal_state
  states_to_windows["INSERT"] = insert_state
  states_to_windows["DELETE"] = delete_state


  print("Number of NORMAL windows: ", len(normal_state))
  print("Number of DELETE windows: ", len(delete_state))
  print("Number of INSERT windows: ", len(insert_state))


  state_to_dist = {}

  for name in states_to_windows:
    state_to_dist[name] = fit_distribution(data=states_to_windows[name],
                                 dist_name=configuration["fit_states_dist"][name])

  # the states of the model.
  # We also need to specify the the probability
  # distribution of the state. this is \pi from the literature
  states = []

  for state in configuration["states"]:
    states.append(State(state_to_dist[state], name=state))


  # add the states to the model
  hmm_model.add_states(states)

  # construct the transition matrix.
  # We create a dense HMM with equal
  # transition probabilities between each state
  # this will be used for initialization when
  # we fit the model. All states have an equal
  # probability to be the starting state or we could

  for i, state in enumerate(configuration["states"]):

    hmm_model.add_transition(hmm_model.start,
                             states[i],
                             configuration["HMM"]["start_prob"][state])

  for i in states:
    for j in states:

      if i == j:
        hmm_model.add_transition(i, j, 0.95)
      else:
        hmm_model.add_transition(i, j, 0.05)

  # finally we need to bake
  hmm_model.bake(verbose=True)
  return hmm_model


def save_windows(windows, configuration, win_interval_length):

  if configuration["save_windows"]:
    import json
    with open(configuration["windows_filename"]+
                  "_"+str(win_interval_length)+".json", 'w') as jsonfile:
      json_str = windows_to_json(windows)
      json.dump(json_str, jsonfile)

def save_hmm(hmm_model, configuration, win_interval_length):

  if configuration["HMM"]["save_model"]:
    json_str = hmm_model.to_json()
    import json
    with open(configuration["HMM"]["save_hmm_filename"]+
              "_"+str(win_interval_length)+".json", 'w') as jsonfile:
      json.dump(json_str, jsonfile)


def zscore_hmm(windows, configuration):

  for win_interval_length in configuration["repeated_sizes"]:

          print("Window interval length: ", win_interval_length)
          n = win_interval_length
          cutoff = (n*configuration["clusterer"]["fpr"]/L)**(1/n)

          print("cutoff used: ", cutoff)

          windows = cluster(data=windows, nclusters=3,
                            method="zscore",
                            **{'n_consecutive_windows': win_interval_length,
                               "cutoff":cutoff})

          # create the HMM
          hmm_model = HiddenMarkovModel(name=configuration["HMM"]["name"],
                                    start=None, end=None)

          # initialize the model
          hmm_model = init_hmm(hmm_model=hmm_model,
                               windows=windows,
                               configuration=configuration)


          flatwindows = flat_windows(windows)

          # fit the model
          hmm_model, history = hmm_model.fit(sequences=flatwindows,
                                         min_iterations=10,
                                         algorithm=configuration["HMM"]["train_solver"],
                                         return_history=True)

          save_hmm(hmm_model=hmm_model,
                   configuration=configuration,
                   win_interval_length=win_interval_length)

def main():

    # load the configuration
    description = "Check the README file for \
      information on how to use the script"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str,
                        default="config_synthetic.json",
                        help='You must specify a json \
                          formatted configuration file')

    args = parser.parse_args()
    config_file = args.config
    configuration = read_configuration_file(config_file)

    windows = create_synthetic_data(configuration=configuration)

    # how many windows in total
    L = len(windows)

    # cluster the windows and assing them state
    if configuration["clusterer"]["name"] =="zscore":
      zscore_hmm(windows=windows, configuration=configuration)

    elif configuration["clusterer"]["name"] =="wmode":
        windows = cluster(data=windows,
                          nclusters=3,
                          method="wmode",
                          **{'normal_rd': configuration["normal_rd"],
                             "delete_rd": configuration["delete_rd"],
                            "insert_rd": configuration["insert_rd"]})
    elif configuration["clusterer"]["name"] =="kmedoids":

        kwargs = {"clusterer":{
                              "name":"kmedoids",
                              "config":{
                              "init_cluster_idx":[0, 10, 15],
                              "metric":"MANHATAN"
                            }
                          }}
        clusterer =  cluster(data=windows,
                          nclusters=3,
                          method="kmedoids",
                          **kwargs)

        from pyclustering.cluster import cluster_visualizer
        from pyclustering.cluster import cluster_visualizer_multidim
        visualizer = cluster_visualizer()
        clusters = clusterer.get_clusters()

        print(clusters)
        clusters_means = []
        cluster_data = defaultdict(list)


        for i in range(len(clusters)):

          window_indexes = clusters[i]

          for w in window_indexes:
            cluster_data[i].extend(windows[i].get_rd_observations())


        for cidx in cluster_data:
          clusters_means.append((cidx, np.mean(cluster_data[cidx])))

        print("Clusters means are: ", clusters_means)



        #visualizer.append_clusters(clusters,
        #                           windows_rd_statistics(windows=windows,
        #                                                 statistic="mean"))
        #visualizer.show()

        #save_windows(windows=windows,
        #           configuration=configuration,
        #           win_interval_length=0)


        """
        # create the HMM
        hmm_model = HiddenMarkovModel(name=configuration["HMM"]["name"],
                                    start=None, end=None)

        # initialize the model
        hmm_model = init_hmm(hmm_model=hmm_model,
                             windows=windows,
                             configuration=configuration)


        flatwindows = flat_windows_from_state(windows=windows,
                                              configuration=configuration,
                                              as_on_seq=False)

        print(flatwindows)

        #flatwindows = flat_windows(windows)
        print("Start training HMM")

        # fit the model
        hmm_model, history = hmm_model.fit(sequences=[flatwindows],
                                           #min_iterations=,
                                           algorithm=configuration["HMM"]["train_solver"],
                                           return_history=True,
                                           verbose=True,
                                           lr_decay=0.6,
                                           callbacks=[HMMCallback(callback=print_logs_callback)],
                                           inertia=0.01)

        print("Done training HMM")


        hmm_model.bake()
        synthetic_data = create_synthetic_data(configuration=configuration,
                                               create_windows=False)

        p_d_given_m = hmm_model.log_probability(sequence=flatwindows)

        print("P(D|M): ", p_d_given_m)

        print(hmm_model.predict_proba(flatwindows))
        viterbi_path=hmm_model.viterbi(flatwindows)

        trans, ems = hmm_model.forward_backward( flatwindows )
        print(trans)

        # plot the model
        plt.figure( figsize=(10,6) )
        hmm_model.plot()
        plt.show()


        save_hmm(hmm_model=hmm_model,
                   configuration=configuration,
                   win_interval_length=0)
        """



if __name__ == '__main__':
  main()

