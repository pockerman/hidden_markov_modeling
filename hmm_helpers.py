"""
Helpers for HMM
"""

import json
from pomegranate import *

from helpers import INFO
from exceptions import Error


class HMMCallback(object):

    def __init__(self, callback):
        self._callback = callback
        self.model = None

    def on_epoch_end(self, logs):
        for key in logs.keys():
            print("{0} {1}: {2}".format(INFO, key, logs[key]))

    def on_training_begin(self):
        pass

    def on_training_end(self, logs):
        pass


def get_window_ids_from_viterbi_path(path, wstate, limit_state):
    indices = []
    for i, item in enumerate(path):

        if item[1].name == wstate:
            indices.append((i, item[1].name))
        elif item[1].name == limit_state:
            indices.append((i, item[1].name))

    return indices


def save_hmm(hmm_model, configuration, win_interval_length):
    json_str = hmm_model.to_json()
    import json
    with open(configuration["HMM"]["save_hmm_filename"] +
              "_" + str(win_interval_length) + ".json", 'w') as jsonfile:
        json.dump(json_str, jsonfile)


def build_hmm(hmm_file):
    with open(hmm_file) as json_file:
        hmm_json_map = json.load(json_file)
        hmm_json_map = json.loads(hmm_json_map)

        hmm = HiddenMarkovModel(name=hmm_json_map["name"],
                                start=None, end=None)

        # reade in the states
        states = hmm_json_map["states"]
        distribution_ties = hmm_json_map.get("distribution ties", None)
        hmm, hmm_states = build_states(hmm=hmm,
                                       states=states,
                                       distribution_ties=distribution_ties)

        hmm.start = hmm_states[hmm_json_map['start_index']]
        hmm.end = hmm_states[hmm_json_map['end_index']]

        # Add all the edges to the model
        for start, end, probability, pseudocount, group in hmm_json_map['edges']:
            hmm.add_transition(hmm_states[start], hmm_states[end], probability,
                               pseudocount, group)

        hmm.bake(verbose=True)
        return hmm


def build_states(hmm, states, distribution_ties):
    state_objs = []
    for state in states:

        state_obj = build_state(state_map=state)

        if state_obj is not None:
            state_objs.append(state_obj)

    if distribution_ties is not None:
        for i, j in distribution_ties:
            # Tie appropriate states together
            states[i].tie(states[j])
        hmm.add_states(state_objs)
    return hmm, state_objs


def build_state(state_map):
    name = state_map["name"]
    weight = state_map["weight"]
    dist_map = state_map["distribution"]

    if dist_map is not None:

        # the state has IndependentComponentsDistribution
        # as a distribution. In this case we have more
        # than one parameters unless we deal with a GMM
        # that wraps the components
        if "name" in dist_map and \
                dist_map["name"] == "IndependentComponentsDistribution":
            parameters = dist_map["parameters"]

            dist_param = parameters[0]
            components = []

            for param in dist_param:

                if param["class"] == "GeneralMixtureModel":

                    # get the dists list for the GMM
                    distributions = param["distributions"]
                    dist_list = []

                    for dist in distributions:
                        distribution = Distribution.from_json(json.dumps(dist))
                        dist_list.append(distribution)

                    weights = param["weights"]
                    gmm = GeneralMixtureModel(dist_list, weights=weights)
                    components.append(gmm)
                elif param["class"] == "Distribution":
                    distribution = Distribution.from_json(json.dumps(param))
                    components.append(distribution)

            # now that we collected the independent components
            # construct the state
            return State(IndependentComponentsDistribution(components),
                         name=name, weight=weight)
        elif "class" in dist_map and \
                dist_map["class"] == "GeneralMixtureModel":

            # get the distributions list for this  GMM
            distributions = dist_map["distributions"]
            dist_list = []

            for dist in distributions:
                distribution = Distribution.from_json(json.dumps(dist))
                dist_list.append(distribution)

            weights = dist_map["weights"]
            gmm = GeneralMixtureModel(dist_list, weights=weights)

            return State(gmm, name=name, weight=weight)
        elif "class" in dist_map and \
                dist_map["class"] == "Distribution":
            distribution = Distribution.from_json(json.dumps(dist_map))
            return State(distribution, name=name, weight=weight)
    else:

        # this means that the state has
        return State(None, name=name, weight=weight)
