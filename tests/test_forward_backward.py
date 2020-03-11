"""
Unit tests for forward-backward algorithms
"""
import unittest
import numpy as np
from forward_backward import forward
from forward_backward import backward

class TestForwardBackward(unittest.TestCase):

    def test_forward(self):
        """
        Test forward algorithm.
        """
        observations = np.array([3, 1, 3])

        # transition matrix
        A = np.array([[0.7, 0.3],
                           [0.4, 0.6]])
        #emission matrix
        B  =  np.array([[0.2, 0.4, 0.4],
                         [0.5, 0.4, 0.1]])

        # initial probability
        initialprob = np.array([0.8, 0.2])

        obs_indices = {1: 0, 2: 1, 3: 2}

        alpha = forward(obs_indices=obs_indices, observations=observations,
                        A=A, B=B, initial_dist=initialprob)

    def test_backward(self):
        """
        Test backward algorithm.
        """
        observations = np.array([3, 1, 3])

        # transition matrix
        A = np.array([[0.7, 0.3],
                      [0.4, 0.6]])
        # emission matrix
        B = np.array([[0.2, 0.4, 0.4],
                      [0.5, 0.4, 0.1]])

        obs_indices = {1: 0, 2: 1, 3: 2}
        beta = backward(obs_indices=obs_indices, dataset=observations, A=A, B=B)


if __name__ == '__main__':
    unittest.main()