import argparse
import pysam
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pomegranate import*

from helpers import read_configuration_file
from helpers import set_up_logger
from bam_helpers import extract_windows
from preprocess_utils import fit_distribution
from exceptions import Error


def main():

    print("Starting analysis")
    description = "Check the README file for information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help='You must specify a json formatted configuration file')
    args = parser.parse_args()

    config_file = args.config
    configuration = read_configuration_file(config_file)

    # configure the logger to use
    set_up_logger(configuration=configuration)
    logging.info("Checking if logger is sane...")

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

    try:

        # TODO: Extractig initial windows is independent
        # we can do so in parallel

        # extract the windows for the WGA treated file
        wga_windows = extract_windows(chromosome=chromosome,
                                      ref_filename=configuration["reference_file"]["name"],
                                      test_filename=configuration["test_file"]["filename"], **args)

        if len(wga_windows) == 0:
            raise Error("WGA windows have not been created")
        else:
            print("\t\tNumber of windows: ", len(wga_windows))

        non_wga_start_idx = configuration["no_wga_file"]["start_idx"]
        non_wga_end_idx = configuration["no_wga_file"]["end_idx"]

        args = {"start_idx": int(non_wga_start_idx),
                "end_idx": (non_wga_end_idx),
                "windowsize": int(windowsize)}

        # exrtact the non-wga windows
        non_wga_windows = extract_windows(chromosome=chromosome,
                                          ref_filename=configuration["reference_file"]["name"],
                                          test_filename=configuration["no_wga_file"]["filename"], **args)

        if len(non_wga_windows) == 0:
            raise Error("Non-WGA windows have not  been created")
        else:
            print("\t\tNumber of windows: ", len(wga_windows))

    except KeyError as e:
        logging.error("Key: {0} does not exit".format(str(e)))
    except Error as e:
        logging.error(str(e))
    except Exception as e:
        logging.error("Unknown exception occured: " + str(e))

    print("Extracted dataset....")
    print("Finished analysis")

    sns.set(color_codes=True)

    # accumulate all RD observations
    wga_rd_observations = []
    non_wga_rd_observations = []

    # now that we have the windows lets do some EDA
    print("Number of windows: ", len(wga_rd_observations))

    # get basic statistics from the windos
    for idx, window in enumerate(wga_windows):
        print("Window id ", idx)
        stats = window.get_rd_stats(statistics="all")
        print("\t mean: ",   stats["mean"])
        print("\t var: ",    stats["var"])
        print("\t median: ", stats["median"])

        rd_data = window.get_rd_observations()
        wga_rd_observations.extend(rd_data)
        sns.distplot(rd_data)
        plt.show()

    mean = np.mean(wga_rd_observations)
    var = np.var(wga_rd_observations)
    median = np.median(wga_rd_observations)

    print("mean: ", mean)
    print("var: ", var)
    print("median: ", median)

    sns.distplot(wga_rd_observations)
    plt.show()

    # get basic statistics from the windos
    for idx, window in enumerate(non_wga_windows):
        print("Window id ", idx)
        stats = window.get_rd_stats(statistics="all")
        print("\t mean: ", stats["mean"])
        print("\t var: ", stats["var"])
        print("\t median: ", stats["median"])

        rd_data = window.get_rd_observations()
        non_wga_rd_observations.extend(rd_data)
        sns.distplot(rd_data)
        plt.show()

    mean = np.mean(non_wga_rd_observations)
    var = np.var(non_wga_rd_observations)
    median = np.median(non_wga_rd_observations)

    print("mean: ", mean)
    print("var: ", var)
    print("median: ", median)

    sns.distplot(non_wga_rd_observations)
    plt.show()

    wga_gaussian_dist = fit_distribution(data=wga_rd_observations)
    wga_gaussian_dist.plot(1000, edgecolor='c', color='c', bins=20)
    plt.show()
    #print("WGA  dist mean: ", wga_gaussian_dist.mu)
    #print("WGA dist var: ", wga_gaussian_dist.cov)

    non_wga_gaussian_dist = fit_distribution(data=non_wga_rd_observations)

    non_wga_gaussian_dist.plot(1000, edgecolor='c', color='c', bins=20)
    plt.show()
    #print("WGA  dist mean: ", non_wga_gaussian_dist.mu)
    #print("WGA dist var: ", non_wga_gaussian_dist.cov)

    # uniform distribution generated
    # from the available data
    #uniform_dist = UniformDistribution.from_samples(rd_observations)
    #uniform_dist.plot(1000, edgecolor='c', color='c', bins=20)
    #plt.show()
    #print("Uniform dist mean: ", uniform_dist.mu)
    #print("Uniform dist var: ", uniform_dist.cov)









if __name__ == '__main__':
    main()