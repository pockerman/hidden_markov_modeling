import argparse
import pysam
from _helpers import read_configuration_file


def main():

    print("Starting analysis")
    description = "Check the README file for information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='config.json',
                        help='You must specify a json formatted configuration file')
    args = parser.parse_args()

    config_file = args.config
    configuration = read_configuration_file(config_file)
    print("Configuration: ", configuration)

    try:
        samfile = pysam.AlignmentFile('m605_verysensitive_trim_sorted.bam',"rb")
    except Exception as e:
        print( str(e))

    print("Finished analysis")

if __name__ == '__main__':
    main()