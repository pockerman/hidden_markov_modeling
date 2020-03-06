"""
Unit tests for HMM class
"""

import unittest
import pandas as pd
import numpy as np
from baumwelch import BaumWelchInput
from baumwelch import BaumWelch

class TestHMM(unittest.TestCase):

    def test_fit(self):
        """
        Test Scenario:
        Expected Output:
        """

        # EqualProbabilities for the initial distribution
        initial_distribution = np.array((0.5, 0.5))
        input = BaumWelchInput(n_itrs=100, show_itrs=True, tol=None, initial_dist=initial_distribution)

        dataset = pd.read_csv("datasets/data_python.csv")
        dataset = dataset["Visible"].values

        # Transition Probabilities
        A = np.ones((2, 2))
        A = A / np.sum(A, axis=1)

        # Emission Probabilities
        B = np.array(((1, 3, 5), (2, 4, 6)))
        B = B / np.sum(B, axis=1).reshape((-1, 1))

        bw  = BaumWelch(input=input)
        bw.solve(dataset=dataset, A=A, B=B)

if __name__ == '__main__':
    unittest.main()
