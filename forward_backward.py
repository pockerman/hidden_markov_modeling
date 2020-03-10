"""
Forward-Backward algos
"""

import numpy as np

def forward(dataset, A, B, initial_dist):
    """
    Forward step.
    See also: http://www.katrinerk.com/courses/python-worksheets/demo-the-forward-backward-algorithm
    :param dataset: Observations
    :param A: the transition probability matrix
    :param B: the observation probability matrix
    :param initial_dist:
    :return:
    """
    # how many observations we have
    n_observations = dataset.shape[0]

    # infer the number of states from the shape
    # of the transition matrix
    n_states = A.shape[0]

    alpha = np.zeros((n_states, n_observations))
    alpha[0, :] = initial_dist * B[:, dataset[0]]

    # the iteration step
    for t in range(1, n_observations):
        for state in range(n_states):
            alpha[t, state] = alpha[t - 1].dot(A[:, state]) * B[state, dataset[t]]

    return alpha