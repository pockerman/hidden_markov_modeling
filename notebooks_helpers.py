"""
Utility helpers used in the stories/ notebooks
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from sklearn import mixture
from sklearn import metrics
from pomegranate import *
from scipy.stats import kde
import csv


def load_data_file(filename):
    """
    Loads a .txt data file into an array
    """
    with open(filename) as file:
        context = file.read()
        size = len(context)
        arraystr= context[1:size-1]
        arraystr = arraystr.split(',')
        region_means = [float(item) for item in arraystr]
        return region_means


def make_data_array(wga_mu, no_wga_mu, gc, use_ratio, use_gc):
    data = []

    if use_ratio and use_gc:
        for no_wga_val, wga_val, gc_val in zip(no_wga_mu, wga_mu, gc):
            data.append([no_wga_val, wga_val, (wga_val + 1) / (no_wga_val + 1), gc_val])
    elif use_ratio:
        for no_wga, wga in zip(no_wga_mu, wga_mu):
            data.append([no_wga, wga, (wga + 1) / (no_wga + 1)])
    elif use_gc:

        for no_wga_val, wga_val, gc_val in zip(no_wga_mu, wga_mu, gc):
            data.append([no_wga_val, wga_val, gc_val])
    else:

        for no_wga, wga in zip(no_wga_mu, wga_mu):
            data.append([no_wga, wga])

    return data


def gmm_clustering(clusters, data, cov_type, tol, max_itrs,
                   n_init, no_wga_mu, wga_mu, xlim, ylim, nbins=80):

    print("Number of clusters ", clusters)
    gmm = mixture.GaussianMixture(n_components=clusters,
                                  covariance_type=cov_type,
                                  tol=tol, max_iter=max_itrs,
                                  n_init=n_init)
    gmm.fit(data)
    print("Converged: ", gmm.converged_)
    print("BIC: ", gmm.bic(data))
    labels = gmm.predict(data)

    print("Len of labels: ", len(labels))

    colors = np.array(['green', 'blue', 'red',
                       'yellow', 'pink', 'orange', 'purple', 'navy',
                       'brown'])

    # add black color for outliers (if any)
    colors = np.append(colors, ["#000000"])
    colors = colors[labels]

    plt.scatter(no_wga_mu, wga_mu, color=colors)
    plt.xlabel("NO-WGA ")
    plt.ylabel("WGA")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()

    # map that holds the association between the
    # component color and component index

    color_comp_assoc = {}
    for label, color in zip(labels, colors):
        if color in color_comp_assoc.keys():
            assert color_comp_assoc[color][0] == label
            color_comp_assoc[color][1] += 1
        else:
            color_comp_assoc[color] = [label, 1]

    return gmm, labels, color_comp_assoc


def is_solid_subseq(subseq, wstate):
    itis = True
    for sitem in subseq:
        if sitem[1].name != wstate:
            itis = False
            break
    return itis


def clean_up_subsequence(subseq, start):
    indices = []

    for item in subseq:
        indices.append((start, item[1].name))
        start += 1
    return indices


def get_state_bwd(path, counter_end, wstate):
    indices = []
    global_end = counter_end
    for item in reversed(path):
        if item[1].name == wstate:
            indices.append((global_end, wstate))
        else:
            break
        global_end -= 1
    if len(indices) != 0:
        return indices[::-1]
    return indices


def get_state_fwd(path, counter_start, wstate):
    indices = []
    global_start = counter_start
    for item in path:

        if item[1].name == wstate:
            indices.append((global_start, item[1].name))
        else:
            break
        global_start += 1
    return indices


def check_items(included, check_on):
    if check_on[-1][0] in included:
        # if the last is included all
        # the rest should be
        return True
    # are not in so put it in
    for item in check_on:
        included.append(item[0])
    return False


def filter_viterbi_path(path, wstate, limit_state, min_subsequence):
    print("Length of path ", len(path))

    indices = []
    index_included = []
    path_counter = 0
    stop = False
    while not stop:

        if path_counter >= len(path):
            break

        # get the item from the viterbi path
        item = path[path_counter]

        # if this is the last item treat it differently
        # as we don't have any forward items
        if path_counter == len(path) - 1:
            # regardless of what this is ignore it
            break

        if item[1].name == limit_state:

            # we don't want to fall off the array
            if min_subsequence + path_counter >= len(path):
                break

            # get the subsequence
            subsequence = path[path_counter:min_subsequence + path_counter]

            # only do work if it is pure
            if is_solid_subseq(subseq=subsequence, wstate=limit_state):

                # bring the subsequence in the form (idx, state)
                subsequence = clean_up_subsequence(subseq=subsequence, start=path_counter)

                # what is after the extracted region
                counter_after = path_counter + min_subsequence

                # check what is after the extracted region
                has_wstate_after = False
                has_limit_state_after = False

                if path[counter_after][1].name == wstate:
                    # this is an island
                    has_wstate_after = True
                elif path[counter_after][1].name == limit_state:
                    has_limit_state_after = True

                # check what is before the extracted
                # region
                has_wstate_before = False

                if path[path_counter - 1][1].name == wstate:
                    has_wstate_before = True

                if has_wstate_before and has_wstate_after:
                    # subsequence is an island get the
                    # surrounding TUFs
                    tuf_before = get_state_bwd(path=path[:path_counter],
                                               counter_end=path_counter - 1,
                                               wstate=wstate)
                    tuf_after = get_state_fwd(path=path[counter_after:],
                                              counter_start=counter_after,
                                              wstate=wstate)
                    if not check_items(included=index_included, check_on=tuf_before):
                        indices.extend(tuf_before)
                    indices.extend(subsequence)
                    indices.extend(tuf_after)

                    # TUF after has been considered
                    # so record it
                    for titem in tuf_after:
                        index_included.append(titem[0])

                    path_counter += min_subsequence + len(tuf_after)

                elif has_wstate_before and has_limit_state_after:
                    # pick up all the limit_wstate as these are part
                    # of the island
                    limit_state_after = get_state_fwd(path=path[counter_after:],
                                                      counter_start=counter_after,
                                                      wstate=limit_state)
                    counter_after += len(limit_state_after)

                    # only if the exactly next one is wstate
                    # then we have an island
                    if path[counter_after][1].name == wstate:

                        tuf_before = get_state_bwd(path=path[:path_counter],
                                                   counter_end=path_counter - 1,
                                                   wstate=wstate)

                        tuf_after = get_state_fwd(path=path[counter_after:],
                                                  counter_start=counter_after,
                                                  wstate=wstate)
                        subsequence.extend(limit_state_after)

                        if not check_items(included=index_included, check_on=tuf_before):
                            indices.extend(tuf_before)

                        indices.extend(subsequence)
                        indices.extend(tuf_after)

                        # TUF after has been considered
                        # so record it
                        for titem in tuf_after:
                            index_included.append(titem[0])

                        path_counter += len(subsequence) + len(tuf_after)
                    else:
                        # jump to counter_after
                        # as did a search and there are no islands
                        path_counter = counter_after
                else:

                    path_counter += 1
            else:
                path_counter += 1
        else:

            # simply increase the counter
            path_counter += 1

    return indices


def get_continuous(tuf_delete_tuf, start_tuf_counter, name):
    has_more = True
    item = None
    counter = start_tuf_counter
    while has_more:

        if counter >= len(tuf_delete_tuf):
            break

        new_item = tuf_delete_tuf[counter]

        if new_item[1] != name:
            has_more = False

        if item is not None:
            if int(new_item[0]) != int(item[0]) + 1:
                has_more = False

        if has_more:
            item = new_item
            counter += 1
    return item, counter


def get_start_end_segment(tuf_delete_tuf, sequence):

    start_tuf_counter = 0
    segments = []
    while True:

        start_item = tuf_delete_tuf[start_tuf_counter]
        start_seq = sequence[start_item[0]]

        end_item, counter = get_continuous(tuf_delete_tuf=tuf_delete_tuf,
                                           start_tuf_counter=start_tuf_counter + 1,
                                           name=start_item[1])
        if end_item is not None:
            assert end_item[1] == start_item[1]
            end_seq = sequence[end_item[0]]

            gap = int(end_seq[1][1]) - int(start_seq[1][0]) + 1
            segments.append((start_item[0], int(start_seq[1][0]),
                             int(end_seq[1][1]), gap,
                             start_item[1]))
        else:
            # this is on its own
            gap = int(start_seq[1][1]) - int(start_seq[1][0])
            segments.append((start_item[0], int(start_seq[1][0]), int(start_seq[1][1]), gap + 1, start_item[1]))

        start_tuf_counter = counter

        if start_tuf_counter >= len(tuf_delete_tuf):
            break

    return segments


def save_segments(segments, chromosome, filename):

    with open(filename, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')

        for segment in segments:
            row = [chromosome]
            row.extend(list(segment))
            writer.writerow(row)


def create_viterbi_path(sequence, hmm_model, filename):
    observations = []
    for i in range(len(sequence)):
        observations.append(sequence[i][0])

    print("Observation length: ", len(observations))

    time_start = time.perf_counter()
    viterbi_path = hmm_model.viterbi(observations)
    time_end = time.perf_counter()
    print("Done. Execution time"
          " {0} secs".format(time_end - time_start))
    print("Log-probability of ML Viterbi path: ", viterbi_path[0])

    # for each item in the sequence
    # cache the index and state predicted
    sequence_viterbi_state = []

    if viterbi_path[1] is not None:
        print("Viterbi path length: ", len(viterbi_path[1]))

        counter = 0
        with open(filename, 'w') as f:
            f.write(str(len(viterbi_path[1]) - 1) + "\n")
            for item in range(len(sequence)):

                if sequence[item][0] == (-999.0, -999.0):
                    counter += 1

                r = (int(sequence[item][1][0]), int(sequence[item][1][1]))
                name = viterbi_path[1][item + 1][1].name
                f.write(str(item) + ":" + str(r) + ":" + str(sequence[item][0]) + ":" + name + "\n")
                sequence_viterbi_state.append((item, viterbi_path[1][item + 1][1].name))

        print("There should be {0} gaps".format(counter))
    else:
        print("Viterbi path is impossible for the given sequence")

    return viterbi_path, observations, sequence_viterbi_state


def plot_state(state_dist, sample_size, min_, max_, n_bins):

     samples = state_dist.sample(n=sample_size)
     xdata = samples[:,[0]]
     ydata = samples[:,[1]]

     x = [item[0] for item in xdata]
     y = [item[0] for item in ydata]

     plt.hist2d(y, x,
                bins=[n_bins, n_bins], cmap='Blues', density=False,
    # cmax=1000,
    # cmin=0,
                alpha=0.99,
                range=((min_, max_), (min_, max_)))

     plt.show()
     state_dist.plot(bins=n_bins)
     plt.show()


def create_states(states_map, means, covariances, plot=True):

    states={}
    for state in states_map:

        print("state name: ", state)
        name = state
        idx = states_map[state]

        # change the order of the means to match the order of the
        # data we will retrieve below in the prediction step
        mu_no_wga = means[idx][0]
        mu_wga = means[idx][1]
        mu = np.array([mu_wga, mu_no_wga])
        print("State means: ", means[idx])
        cov = covariances[idx]
        print("State covariance: ", cov)
        cov = np.array([[cov[1], 0.0], [0.0, cov[0]]])

        if state == "Duplication":
            mu[1] = 45.0

        state_dist = MultivariateGaussianDistribution(means=mu,
                                                      covariance=cov)

        if plot:

            if state == "Duplication":
                min_ = 0.0
                max_ = 70.0
            elif state == "Normal-I":
                min_ = 0.0
                max_ = 70.0
            elif state == "Normal-II":
                min_ = 0.0
                max_ = 70.0
            elif state == 'Deletion':
                min_ = 0.0
                max_ = 70.0

            plot_state(state_dist=state_dist,
                       sample_size=10000, n_bins=30,
                       min_=min_, max_=max_)

        states[name] = State(state_dist, name=name)
    return states


def create_tuf_state(tuf_means, tuf_cov):

    print("TUF means: ", tuf_means)
    print("TUF covariances: ", tuf_cov)

    tuf_mu_no_wga = tuf_means[0]
    tuf_mu_wga = tuf_means[1]
    tuf_mu = np.array([tuf_mu_wga, tuf_mu_no_wga])
    tuf_cov = tuf_cov
    tuf_cov = np.array([[tuf_cov[1], 0.0], [0.0, tuf_cov[0]]])

    tuf_dist = MultivariateGaussianDistribution(means=tuf_mu,
                                                covariance=tuf_cov)

    n_bins = 100
    samples = tuf_dist.sample(n=10000)
    xdata = samples[:, [0]]
    ydata = samples[:, [1]]

    x = [item[0] for item in xdata]
    y = [item[0] for item in ydata]
    # print(samples)
    plt.hist2d(y, x,
               bins=[n_bins, n_bins], cmap='Blues', density=False,
               cmax=1000,
               cmin=0,
               alpha=0.99,
               range=((0.0, 70.), (0.0, 70.0)))

    plt.show()
    tuf_dist.plot(bins=n_bins)
    plt.show()

    # only change the means
    tuf_mu = np.array([tuf_mu_wga, 40.0])

    # also add the dist that is not modeled by the data
    tuf_dist_double = MultivariateGaussianDistribution(means=tuf_mu,
                                                       covariance=tuf_cov)

    tuf_mixture = GeneralMixtureModel([tuf_dist, tuf_dist_double], weights=[0.66, 0.34])
    return State(tuf_mixture, name='TUF')


def get_states_counter(states_counter, observations, sequence_viterbi_state):

    # collect counters for error statistics
    for i, obs in enumerate(observations):

        # what dows the Viterbi path says
        viterbi_state = sequence_viterbi_state[i][1]
        if obs == (-999.0, -999.0):

            # if we predict that increase the success
            # otherwise increase errors
            if viterbi_state == 'GAP_STATE':
                states_counter[viterbi_state] += 1
            else:
                states_counter["GAP_STATE_INCORRECT"] += 1
        else:
            states_counter[viterbi_state] += 1
    return states_counter


def plot_hmm_states_to_labels(hmm_states_to_labels, observations, sequence_viterbi_state, no_wga_obs, wga_obs):
    hmm_labels = []

    # collect the labels as these are predicted by the HMM
    for i, obs in enumerate(observations):

        # do not account for GAPs
        if obs != (-999.0, -999.0):
            viterbi_state = sequence_viterbi_state[i][1]
            hmm_labels.append(hmm_states_to_labels[viterbi_state])

    colors = np.array(['green', 'blue', 'red', 'yellow', 'pink', 'black'])
    colors = colors[hmm_labels]

    plt.scatter(no_wga_obs, wga_obs, color=colors)

    plt.xlabel("NO-WGA ")
    plt.ylabel("WGA")
    plt.xlim(0.0, 80.0)
    plt.ylim(0.0, 80.0)
    plt.show()

    color_comp_assoc_hmm = {}
    for label, color in zip(hmm_labels, colors):
        if color in color_comp_assoc_hmm.keys():
            assert color_comp_assoc_hmm[color][0] == label
            color_comp_assoc_hmm[color][1] += 1
        else:
            color_comp_assoc_hmm[color] = [label, 1]
    return color_comp_assoc_hmm, hmm_states_to_labels, hmm_labels


def plot_hmm_label_state(hmm_states_to_labels, hmm_labels, no_wga_obs, wga_obs, nbins):

    for label in hmm_states_to_labels:

        print("State: ", label)
        label_idx = hmm_states_to_labels[label]

        state_labels = []

        state_no_wga_obs = []
        state_wga_obs = []
        for i, item in enumerate(hmm_labels):
            if item == label_idx:
                state_labels.append(label_idx)
                state_no_wga_obs.append(no_wga_obs[i])
                state_wga_obs.append(wga_obs[i])

        colors = np.array(['green', 'yellow', 'blue', 'red', 'pink', 'purple'])

        if len(state_no_wga_obs) != 0:

            print("Length: {0}".format(len(state_no_wga_obs)))

            colors = colors[state_labels]

            # plot the observations
            plt.scatter(state_no_wga_obs, state_wga_obs, color=colors)

            kernel = kde.gaussian_kde(np.vstack([state_no_wga_obs, state_wga_obs]))
            xi, yi = np.mgrid[min([0.0]):max([70.0]):nbins * 1j,
                     min([0.0]):max([70.0]):nbins * 1j]

            zi = kernel(np.vstack([xi.flatten(), yi.flatten()]))
            plt.contour(xi, yi, zi.reshape(xi.shape), 24)

            plt.xlabel("NO-WGA ")
            plt.ylabel("WGA")

            if label == "Duplication":
                plt.xlim(0.0, 70.0)
                plt.ylim(0.0, 70.0)
            if label == "Normal-I":
                plt.xlim(0.0, 70.0)
                plt.ylim(0.0, 70.0)
            if label == "Normal-II":
                plt.xlim(0.0, 80.0)
                plt.ylim(0.0, 80.0)
            if label == "Deletion":
                plt.xlim(0.0, 70.0)
                plt.ylim(0.0, 70.0)
            if label == "TUF":
                plt.xlim(0.0, 70.0)
                plt.ylim(0.0, 70.0)

            plt.show()
        else:
            print("For state: {0} could not plot empty observations".format(label))

