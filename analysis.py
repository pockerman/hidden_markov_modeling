import argparse
import pysam
from helpers import read_configuration_file
from bam_helpers import extract_windows
from exceptions import Error
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

    with pysam.AlignmentFile(configuration["reference_file"]["name"],"rb") as ref_file:

        print("\t Reference file")
        print("\t", ref_file.filename)
        print("\t", ref_file.description)
        print("\n")

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
                windows = extract_windows(chromosome=chromosome, ref_file=ref_file,
                                      test_file=test_file, **args)

                if(len(windows) == 0):
                    raise Error("No windows have been created")
                else:
                    print("\t\tNumber of windows: ", len(windows))

                print("Finished analysis")

                """
                # apply preprocessing for the windows
                windows = preprocess(windows=windows)
    
                # specify the HMM model
                hmm = HMM(start_transitions=configuration["HMM"]["initial_transitions_p"])
                hmm.fit(dataset=windows, solver=configuration["HMM"]["train_solver"])
                """
            except KeyError as e:
                print("Key: " + str(e) + " does not exist")
            except Error as e:
                print("Error occured: " + str(e))
            except Exception as e:
                print("Unknown exception occured: " + str(e))





if __name__ == '__main__':
    main()