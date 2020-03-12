"""
Forward-Backward algos
"""

import numpy as np

def forward(obs_indices, observations, A, B, initial_dist):

    """
    Forward step.
    See also: http://www.katrinerk.com/courses/python-worksheets/demo-the-forward-backward-algorithm

    :param obs_indices: A dictionary that maps observations to indices
    :param observations: Observations
    :param A: the transition probability matrix
    :param B: the observation probability matrix
    :param initial_dist:
    :return: numpy matrix of size (N_OBSERVATIONS, N_STATES)
    """
    # how many observations we have
    n_observations = observations.shape[0]

    # infer the number of states from the shape
    # of the transition matrix
    n_states = A.shape[0]

    alpha = np.zeros((n_observations, n_states))

    obs_index = obs_indices[observations[0]]
    alpha[0, :] = initial_dist * B[:, obs_index]

    # the iteration step
    for t in range(1, n_observations):
        for state in range(n_states):
            alpha[t, state] = alpha[t - 1].dot(A[:, state]) * B[state, obs_indices[observations[t]]]

    return alpha


def backward(obs_indices, dataset, A, B):
    """
    Backward step

    beta_j(t) = P(O_{t+1}, ..., O_T | q_t = j, HMM)

    :param obs_indices: A dictionary that maps observations to indices
    :param dataset: Observations
    :param A: transition probability matrix
    :param B: the observation probability matrix
    :return: numpy matrix of size (N_OBSERVATIONS, N_STATES)
    """

    n_states = A.shape[0]
    n_observations = dataset.shape[0]

    beta = np.zeros((n_observations, n_states))

    # setting beta(T) = 1
    for state in range(n_states):
        beta[n_observations - 1][state] = 1

    # Loop in backward way from T-1 to
    # Due to python indexing the actual loop will be T-2 to 0
    for t in range(dataset.shape[0] - 2, -1, -1):
        for j in range(A.shape[0]):
            beta[t, j] = (beta[t + 1] * B[:, obs_indices[dataset[t + 1]]]).dot(A[j, :])

    return beta
