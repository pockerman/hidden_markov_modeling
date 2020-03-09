import argparse
import pysam
from helpers import read_configuration_file
from bam_helpers import extract_windows
from preprocess_utils import  preprocess
from hmm import HMM

def main():

    print("Starting analysis")
    description = "Check the README file for information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help='You must specify a json formatted configuration file')
    args = parser.parse_args()

    config_file = args.config
    configuration = read_configuration_file(config_file)

    print("\tConfiguration: ", configuration)

    try:

        # read the refernce  file
        ref_file = pysam.AlignmentFile(configuration["reference_file"],"rb")

        # read the test file
        test_file = pysam.AlignmentFile(configuration["reference_file"],"rb")

        """
        # extract the windows
        windows = extract_windows(chromosome=configuration["chromosome"], ref_file=ref_file,
                                  start_test=test_file, **{"start_test": configuration["start_test"],
                                                           "end_test": configuration["end_test"]})
        # apply preprocessing for the windows
        windows = preprocess(windows=windows)

        # specify the HMM model
        hmm = HMM(start_transitions=configuration["HMM"]["initial_transitions_p"])
        hmm.fit(dataset=windows, solver=configuration["HMM"]["train_solver"])
        """

    except Exception as e:
        print( str(e))

    print("Finished analysis")

if __name__ == '__main__':
    main()