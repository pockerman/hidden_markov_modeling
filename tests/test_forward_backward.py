"""
Unit tests for forward-backward algorithms
"""
import unittest
import numpy as np
from forward_backward import forward

class TestForwardBackward(unittest.TestCase):

    def test_forward(self):
        """
        Test forward algorith.
        :return:
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

        alpha = forward(dataset=observations, A=A, B=B, initial_dist=initialprob)


if __name__ == '__main__':
    unittest.main()