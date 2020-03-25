from pomegranate import*
import matplotlib.pyplot as plt
import argparse

from helpers import read_configuration_file
from helpers import flat_windows
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
            "windowsize": int(windowsize)}

    # extract the windows for the WGA treated file
    wga_windows = extract_windows(chromosome=chromosome,
                                      ref_filename=configuration["reference_file"]["name"],
                                      test_filename=configuration["test_file"]["filename"], **args)

    if len(wga_windows) == 0:
        raise Error("WGA windows have not been created")
    else:
        print("\t\tNumber of windows: ", len(wga_windows))

    return wga_windows


def main():

    # load the configuration
    description = "Check the README file for information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help='You must specify a json formatted configuration file')
    args = parser.parse_args()

    config_file = args.config
    configuration = read_configuration_file(config_file)

    wga_windows = load_data(configuration=configuration)

    # create the HMM
    model = HiddenMarkovModel(name="HMMModel")

    # the states of the model
    insert = State(UniformDistribution(1, 5), name="INSERT")
    delete = State(UniformDistribution(1, 5), name="DELETE")
    normal = State(UniformDistribution(1, 5), name="NORMAL")

    states = [insert, delete,  normal]

    # add the states to the model
    model.add_states(insert, delete,  normal)

    # construct the transition matrix
    # this will be used for initialization when
    # we fit the model. All states have an equal
    # probability to be the starting state or we could
    # initialize using UniformDistribution().sample()
    model.add_transition(model.start, insert, 1.0/len(states))
    model.add_transition(model.start, delete, 1.0/len(states))
    model.add_transition(model.start, normal, 1.0/len(states))

    # create a dense HMM with equal
    # transition probabilities between each state
    for i in states:
        for j in states:
            model.add_transition(i, j, 0.5)

    # finally we need to bake
    model.bake()

    # fit the model
    model, history = model.fit(sequences=flat_windows(wga_windows),
                        algorithm='baum-welch', return_history=True)

    print("History epoch: ", history)
    #print("History epoch: ", history.learning_rate)

    # plot the model
    plt.figure( figsize=(10,6) )
    model.plot()


if __name__ == '__main__':
    main()