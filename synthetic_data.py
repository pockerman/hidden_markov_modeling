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
#import seaborn as sns
import numpy as np
#import matplotlib.pyplot as plt
#from collections import defaultdict
#from scipy import stats

from helpers import read_configuration_file
from helpers import Window
#from helpers import WindowState
from helpers import Observation
from helpers import add_window_observation
from helpers import DUMMY_ID
#from helpers import windows_to_json
from helpers import flat_windows
from helpers import flat_windows_from_state
from helpers import HMMCallback
from helpers import print_logs_callback
#from helpers import windows_rd_statistics
from helpers import save_hmm

#from helpers import flat_windows_rd_from_indexes
from cluster import Cluster
#from cluster import clusters_statistics
from hypothesis_testing import SignificanceTestLabeler
from exceptions import Error
from bam_helpers import DUMMY_BASE
from preprocess_utils import build_clusterer
from preprocess_utils import fit_distribution


def create_synthetic_data(configuration, create_windows=True):

  bases = ['A', 'T', 'C', 'G']
  flip = ["YES", "NO"]
  states = {"NORMAL": configuration["normal_rd"],
            "DELETE": configuration["delete_rd"],
            "INSERT": configuration["insert_rd"]}

  seq_size = configuration["sequence_size"]

  synthetic_data = []
  state_name = np.random.choice(list(states.keys()))

  counter = 0
  for i in range(seq_size):

    if np.random.choice(flip) == 'YES' and \
      counter > configuration["minimum_occurences"][state_name]:
      state_name = np.random.choice(list(states.keys()))
      synthetic_data.append(states[state_name])
      counter = 0
      counter += 1
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


def create_clusters(windows, configuration):

  kwargs = {"clusterer":{ "name":configuration["clusterer"],
                          "config":configuration["clusterer"]["config"]}}


  clusterer =  build_clusterer(data=windows,
                               nclusters=len(configuration["states"]),
                               method="kmedoids",
                               **kwargs)

  clusters_indexes = clusterer.get_clusters()
  clusters = []

  for i in range(len(clusters_indexes)):
    clusters.append(Cluster(id_ = i, indexes=clusters_indexes[i]))

  labeler = SignificanceTestLabeler(clusters=clusters,
                                    windows=windows)

  labeled_clusters = labeler.label(test_config=configuration["labeler"])

  # update windows states
  for state in labeled_clusters:
    cluster = labeled_clusters[state]
    indexes = cluster.indexes

    for idx in indexes:
      windows[idx].set_state(cluster.state)

  return labeled_clusters

def init_hmm(clusters, windows, configuration):

  # create the HMM
  hmm_model = HiddenMarkovModel(name=configuration["HMM"]["name"],
                                start=None, end=None)


  state_to_dist = {}
  states = []
  for cluster in clusters:
    state_to_dist[cluster.state.name] = \
    fit_distribution(data=cluster.get_data_from_windows(windows=windows),
                     dist_name=configuration["fit_states_dist"][cluster.state.name])
    states.append(State(state_to_dist[cluster.state.name], name=cluster.state.name))


  # add the states to the model
  hmm_model.add_states(states)

  # construct the transition matrix.
  # We create a dense HMM with equal
  # transition probabilities between each state
  # this will be used for initialization when
  # we fit the model. All states have an equal
  # probability to be the starting state or we could

  for i, cluster in enumerate(clusters):
    hmm_model.add_transition(hmm_model.start,
                             states[i],
                              configuration["HMM"]["start_prob"][cluster.state.name])

  for i in states:
    for j in states:

      if i == j:
        hmm_model.add_transition(i, j, 0.95)
      else:
        hmm_model.add_transition(i, j, 0.05)

  # finally we need to bake
  hmm_model.bake(verbose=True)
  return hmm_model


def hmm_train(clusters, windows, configuration):

  # initialize the model
  hmm_model = init_hmm(clusters=clusters,
                       windows=windows,
                       configuration=configuration)


  flatwindows = flat_windows_from_state(windows=windows,
                                        configuration=configuration,
                                        as_on_seq=False)

  #flatwindows = flat_windows(windows=windows)

  print("Flatwindows are: ", flatwindows)

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

    clusters = create_clusters(windows=windows, configuration=configuration)

    print("Number of clusters used: {0}".format(len(clusters)))

    for cluster in clusters:
      print("State modelled by cluster {0} is {1}".format(clusters[cluster].cidx,
                                                          clusters[cluster].state.name))

    hmm_train(clusters=clusters.values(), windows=windows, configuration=configuration)


if __name__ == '__main__':
  main()

