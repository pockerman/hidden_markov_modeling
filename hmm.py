"""
Hidden Markov Model general class
"""


class HMM(object):

    def __init__(self, start_prob, start_transitions):
        """
        Build an HMM instance by supplying the starting probability
        vector and the starting transitions matrix
        :param start_prob:
        :param start_transitions:
        """
        self._start_prob = start_prob
        self._start_transitions = start_transitions
        self._model_params = None

    def fit(self, dataset, solver):
        """
        Given the observed data estimate the optimal
        model parameters
        :param dataset:
        :return:
        """
        pass