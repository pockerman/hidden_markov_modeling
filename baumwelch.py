"""
Baum-Welch algorithm for HMM
fitting
"""

import numpy as np
import collections

BaumWelchInput = collections.namedtuple("BaumWelchInput", ["n_itrs", "tol",
                                                           "show_itrs", "initial_dist"])
BaumWelchOutput = collections.namedtuple("BaumWelchOutput", ["n_itrs", "tol"])


class BaumWelch:

    def __init__(self, input):
        self._input = input

    def solve(self, dataset, A, B):
        """
        Solve for the given transition matrix A and
         emission matrix B on the given data set
        """
        T = len(dataset)
        M = A.shape[0]

        for itr in range(self._input.n_itrs):

            if self._input.show_itrs:
                print("Iteration: ", itr)

            alpha = self._predict(dataset=dataset, A=A, B=B)
            beta = self._update(dataset=dataset, A=A, B=B)

            xi = np.zeros((M, M, T-1))

            for t in range(T-1):
                denom = np.dot(np.dot(alpha[t,:].T, A) * B[:, dataset[t + 1]].T, beta[t + 1, :])

                for i in range(M):
                    numerator = alpha[t, i] * A[i, :] * B[:, dataset[t + 1]].T * beta[t + 1, :].T
                    xi[i, :, t] = numerator / denom

            # calculate the probability being at state S_i
            # at time t (Rabiner eq 38
            gamma = np.sum(xi, axis=1)

            # recalculate A matrix. See Rabiner eq 40b
            A = np.sum(xi, 2) / np.sum(gamma, axis=1).reshape((-1, 1))

            # Add additional T'th element in gamma
            gamma = np.hstack((gamma, np.sum(xi[:, :, T - 2], axis=0).reshape((-1, 1))))

            # Recalculate B. See Rabiner eq 40c
            K = B.shape[1]
            denom = np.sum(gamma, axis=1)
            for l in range(K):
                B[:, l] = np.sum(gamma[:, dataset == l], axis=1)

            B = np.divide(B, denom.reshape((-1, 1)))

        return A, B

    def _predict(self, dataset, A, B):

        alpha = np.zeros((dataset.shape[0], A.shape[0]))
        alpha[0, :] = self._input.initial_dist * B[:, dataset[0]]

        for t in range(1, dataset.shape[0]):
            for j in range(A.shape[0]):
                alpha[t, j] = alpha[t - 1].dot(A[:, j]) * B[j, dataset[t]]
        return alpha

    def _update(self, dataset, A, B):

        beta = np.zeros((dataset.shape[0], A.shape[0]))

        # setting beta(T) = 1
        beta[dataset.shape[0] - 1] = np.ones((A.shape[0]))

        # Loop in backward way from T-1 to
        # Due to python indexing the actual loop will be T-2 to 0
        for t in range(dataset.shape[0] - 2, -1, -1):
            for j in range(A.shape[0]):
                beta[t, j] = (beta[t + 1] * B[:, dataset[t + 1]]).dot(A[j, :])

        return beta