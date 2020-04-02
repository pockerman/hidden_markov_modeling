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

from helpers import read_configuration_file
from helpers import Window
from helpers import WindowState
from helpers import Observation
from helpers import add_window_observation
from helpers import DUMMY_ID
from helpers import windows_to_json
from helpers import flat_windows
from exceptions import Error
from bam_helpers import DUMMY_BASE
from preprocess_utils import cluster
from preprocess_utils import fit_distribution


def create_synthetic_data(configuration):

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

    if np.random.choice(flip) == 'YES' and counter >= 10:
      state_name = np.random.choice(list(states.keys()))
      synthetic_data.append(states[state_name])
      counter += 1
    elif counter == configuration["swap_state_freq"]:

      state_name = np.random.choice(list(states.keys()))
      synthetic_data.append(states[state_name])
      counter += 1

    else:
      synthetic_data.append(states[state_name])
      counter += 1


  if len(synthetic_data) != seq_size:
    raise Error("Invalid size of synthetic \
                data {0} not equal to {1}".format(len(synthetic_data), seq_size))


  print(synthetic_data)

  # let's create the windows
  # the returned list of windows
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

    L = len(windows)

    for win_interval_length in configuration["repeated_sizes"]:

      print("Window interval length: ", win_interval_length)

      n = win_interval_length
      cutoff = (n*configuration["clusterer"]["fpr"]/L)**(1/n)

      print("cutoff used: ", cutoff)

      # cluster the windows and assing them state
      windows = cluster(data=windows, nclusters=3,
                        method="zscore",
                        **{'n_consecutive_windows': win_interval_length,
                           "cutoff":cutoff})


      print("Number of windows: ", len(windows))

      for window in windows:
        print("Window id: ", window.get_id())
        print("Window state: ", window.get_state().name)
        print("Window length: ", len(window))

      # if we want to save the windows then do so
      if configuration["save_windows"]:
        import json
        with open(configuration["windows_filename"]+
                  "_"+str(win_interval_length)+".json", 'w') as jsonfile:
          json_str = windows_to_json(windows)
          json.dump(json_str, jsonfile)

      # collect the windows with the same state
      normal_state = []
      delete_state = []
      insert_state = []

      for window in windows:
        if window.get_state() == WindowState.NORMAL:
          normal_state.extend(window.get_rd_observations())
        elif window.get_state() == WindowState.DELETE:
          delete_state.extend(window.get_rd_observations())
        elif window.get_state() == WindowState.INSERT:
          insert_state.extend(window.get_rd_observations())


      normal_dist = fit_distribution(data=normal_state,
                                     dist_name=configuration["fit_states_dist"]["NORMAL"])
      delete_dist = fit_distribution(data=delete_state,
                                     dist_name=configuration["fit_states_dist"]["DELETE"])
      insert_dist = fit_distribution(data=insert_state,
                                     dist_name=configuration["fit_states_dist"]["INSERT"])

      # create the HMM
      model = HiddenMarkovModel(name=configuration["HMM"]["name"],
                                start=None, end=None)

      # the states of the model.
      # We also need to specify the the probability
      # distribution of the state. this is \pi from the literature
      insert = State(insert_dist, name="INSERT")
      delete = State(delete_dist, name="DELETE")
      normal = State(normal_dist, name="NORMAL")

      states = [insert, delete,  normal]

      # add the states to the model
      model.add_states(insert, delete,  normal)

      # construct the transition matrix.
      # We create a dense HMM with equal
      # transition probabilities between each state
      # this will be used for initialization when
      # we fit the model. All states have an equal
      # probability to be the starting state or we could
      # initialize using UniformDistribution().sample()
      model.add_transition(model.start, insert, 1.0/len(states))
      model.add_transition(model.start, delete, 1.0/len(states))
      model.add_transition(model.start, normal, 1.0/len(states))

      for i in states:
          for j in states:
              model.add_transition(i, j, 0.5)

      # finally we need to bake
      model.bake()

      flatwindows = flat_windows(windows)

      # fit the model
      model, history = model.fit(sequences=flatwindows,
                                 min_iterations=10,
                                 algorithm=configuration["HMM"]["train_solver"],
                                 return_history=True)

      # save the model
      if configuration["HMM"]["save_model"]:
        json_str = model.to_json()
        import json
        with open(configuration["HMM"]["save_hmm_filename"]+
                  "_"+str(win_interval_length)+".json", 'w') as jsonfile:
          json.dump(json_str, jsonfile)



if __name__ == '__main__':
  main()

