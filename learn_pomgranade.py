from pomegranate import*
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np

from helpers import read_configuration_file
from helpers import flat_windows
from helpers import windows_to_json
from preprocess_utils import cluster
from preprocess_utils import fit_distribution
from bam_helpers import extract_windows
from exceptions import Error

def load_data(configuration):
    """
    Load the data
    """

    wga_start_idx = configuration["reference_file"]["start_idx"]
    wga_end_idx = configuration["reference_file"]["end_idx"]
    windowsize = configuration["window_size"]
    chromosome = configuration["chromosome"]

    print("\t\tStart index used: ", wga_start_idx)
    print("\t\tEnd index used: ", wga_end_idx)
    print("\t\tWindow size: ", windowsize)
    print("\t\tChromosome: ", chromosome)

    args = {"start_idx": int(wga_start_idx),
            "end_idx": (wga_end_idx),
            "windowsize": int(windowsize),
            "quality_theshold": configuration["quality_theshold"],
            "fill_missing_window_data": configuration.get("fill_missing_window_data", False),
            "fill_missing_window_data_factor": configuration.get("fill_missing_window_data_factor", 0)
          }

    ref_filename = configuration["reference_file"]["name"]
    test_filename = configuration["test_file"]["filename"]

    # extract the windows for the WGA treated file
    wga_windows = extract_windows(chromosome=chromosome,
                                  ref_filename=ref_filename,
                                  test_filename=test_filename,
                                  **args)

    if len(wga_windows) == 0:
        raise Error("WGA windows have not been created")
    else:
        print("\t\tNumber of windows: ", len(wga_windows))

    return wga_windows


def build_hmm():
  pass


def main():

    # load the configuration
    description = "Check the README file for information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help='You must specify a json formatted configuration file')
    args = parser.parse_args()

    config_file = args.config
    configuration = read_configuration_file(config_file)

    # load the windowed data
    wga_windows = load_data(configuration=configuration)

    # if we want to save the windows then do so
    if configuration["save_windows"]:
      import json
      with open(configuration["windows_filename"], 'w') as jsonfile:
        json_str = windows_to_json(wga_windows)
        json.dump(json_str, jsonfile)


    windows = flat_windows(wga_windows)

    # use clustering to cluster similar windows together
    # then use the resulting clusters to fit probability
    # distributions to use for initializing the emission probabilities
    # for each state
    clusterer = cluster(data=windows, nclusters=3, method="kmeans")

    # visualize the clusters
    clusters = clusterer.labels_

    cluster0 = []
    cluster1 = []
    cluster2 = []

    for c in range(len(clusters)):

      if clusters[c] == 0:
        cluster0.extend(windows[c])
      elif clusters[c] == 1:
        cluster1.extend(windows[c])
      elif clusters[c] == 2:
        cluster2.extend(windows[c])
      else:
        raise ValueError("Unknown cluster id")

    dist_clust_0 = fit_distribution(data=cluster0)
    dist_clust_1 = fit_distribution(data=cluster1)
    dist_clust_2 = fit_distribution(data=cluster2)

    # create the HMM
    model = HiddenMarkovModel(name=configuration["HMM"]["name"],
                              start=None, end=None)

    # the states of the model.
    # We also need to specify the the probability
    # distribution of the state. this is \pi from the literature
    insert = State(dist_clust_0, name="INSERT")
    delete = State(dist_clust_1, name="DELETE")
    normal = State(dist_clust_2, name="NORMAL")

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

    # fit the model
    model, history = model.fit(sequences=windows,
                               min_iterations=10,
                               algorithm=configuration["HMM"]["train_solver"],
                               return_history=True)

    # save the model
    if configuration["HMM"]["save_model"]:
      json_str = model.to_json()
      import json
      with open(configuration["HMM"]["save_hmm_filename"], 'w') as jsonfile:
        json.dump(json_str, jsonfile)

    #model.bake()

    #print("Model is: ", model)

    #print("History log: ", history.log)
    #print("History epoch: ", history.learning_rate)

    # plot the model
    #plt.figure( figsize=(10,6) )
    #model.plot()


if __name__ == '__main__':
    main()