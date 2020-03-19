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
    logging.info("Checking logger...")

    with pysam.FastaFile(configuration["reference_file"]["name"]) as ref_file:

        print("\t Reference file: ", ref_file.filename)

        with pysam.AlignmentFile(configuration["test_file"]["filename"],"rb") as test_file:

            print("\t Test file")
            print("\t", test_file.filename)
            print("\t", test_file.description)

            print("=============================\n")

            try:

                print("\t Extracting windows")
                start_idx= configuration["reference_file"]["start_idx"]
                end_idx = configuration["reference_file"]["end_idx"]
                windowsize = configuration["window_size"]
                chromosome = configuration["chromosome"]

                print("\t\tStart index used: ", start_idx)
                print("\t\tEnd index used: ",   end_idx)
                print("\t\tWindow size: ", windowsize)
                print("\t\tChromosome: ", chromosome)

                args = {"start_idx":int(start_idx),
                        "end_idx":(end_idx),
                        "windowsize":int(windowsize)}

                # extract the windows
                windows = extract_windows(chromosome=chromosome,
                                          ref_file=ref_file,
                                          test_file=test_file, **args)

                if len(windows) == 0:
                    raise Error("No windows have been created")
                else:
                    print("\t\tNumber of windows: ", len(windows))

            except KeyError as e:
                logging.error("Key: {0} does not exit".format(str(e)))
            except Error as e:
                logging.error(str(e))
            except Exception as e:
                logging.error("Unknown exception occured: " + str(e))

    print("Extracted dataset....")
    print("Finished analysis")

    """
    sns.set(color_codes=True)

    # accumulate all RD observations
    rd_observations = []

    # now that we have the windows lets do some EDA
    print("Number of windows: ", len(windows))

    # get basic statistics from the windos
    for idx, window in enumerate(windows):
        print("Window id ", idx)
        stats = window.get_rd_stats(statistics="all")
        print("\t mean: ",   stats["mean"])
        print("\t var: ",    stats["var"])
        print("\t median: ", stats["median"])

        rd_data = window.get_rd_observations()
        rd_observations.extend(rd_data)
        sns.distplot(rd_data)
        plt.show()

    mean = np.mean(rd_observations)
    var = np.var(rd_observations)
    median = np.median(rd_observations)

    print("mean: ", mean)
    print("var: ", var)
    print("median: ", median)

    sns.distplot(rd_observations)
    plt.show()

    # uniform distribution generated
    # from the available data
    uniform_dist = UniformDistribution.from_samples(rd_observations)
    uniform_dist.plot(1000, edgecolor='c', color='c', bins=20)
    plt.show()
    #print("Uniform dist mean: ", uniform_dist.mu)
    #print("Uniform dist var: ", uniform_dist.cov)
    """








if __name__ == '__main__':
    main()