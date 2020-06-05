# Test on the Viterbi algorithm
import time
import sys

sys.path.append("../")
import json
import matplotlib.pyplot as plt
import seaborn as sns


from helpers import read_configuration_file
from train import main as train_main
from train import load_regions
from hmm_helpers import build_hmm
from helpers import WindowType


def main():

  # load the configuration
  configuration=read_configuration_file("../config.json")
  configuration["HMM"]["train"]=False


  configuration["regions_files"]=["/home/a/ag568/region_0_MANHATAN_3_MEAN_RATIO.txt"]
  clusters = {
    "cluster_0":{"filename":"EMPTY", "state":"EMPTY", "distributions":{} },
    "cluster_1":{"filename":"EMPTY", "state":"EMPTY", "distributions":{}},
    "cluster_2":{"filename":"EMPTY", "state":"EMPTY", "distributions":{}}
    }

  clusters["cluster_0"]["filename"]="/home/a/ag568/cluster_0_MANHATAN_3_MEAN_RATIO.txt"
  clusters["cluster_0"]["state"]="state_0"

  clusters["cluster_1"]["filename"]="/home/a/ag568/cluster_1_MANHATAN_3_MEAN_RATIO.txt"
  clusters["cluster_1"]["state"]="state_1"

  clusters["cluster_2"]["filename"]="/home/a/ag568/cluster_2_MANHATAN_3_MEAN_RATIO.txt"
  clusters["cluster_2"]["state"]="state_2"


  configuration["clusters"] = clusters

  hmm_config = configuration["HMM"]

  start_prob = 0.25
  hmm_config["states"]= {"state_0":{ "start_prob":start_prob},
                         "state_1":{ "start_prob":start_prob},
                         "state_2":{ "start_prob":start_prob},
                         "gap_state": {"start_prob":start_prob}
                         }

  self_trans = 0.95
  off_trans = 0.016
  hmm_config["transitions"]={
      "state_0-state_0":self_trans,
      "state_0-state_1":off_trans,
	  "state_0-state_2":off_trans,
      "state_0-gap_state":off_trans,
      "state_1-state_0":off_trans,
      "state_1-state_1":self_trans,
      "state_1-state_2":off_trans,
      "state_1-gap_state":off_trans,
      "state_2-state_0":off_trans,
      "state_2-state_1":off_trans,
      "state_2-state_2":self_trans,
      "state_2-gap_state":off_trans,
      "gap_state-state_0":off_trans,
      "gap_state-state_1":off_trans,
      "gap_state-state_2":off_trans,
      "gap_state-gap_state":self_trans,

    }

  # now we can train
  hmm, regions = train_main(configuration=configuration)
  print("Number of regions: ", len(regions))

  for region in regions:
      print("Number of gaps in region: {0} is {1}".format(region.ridx, region.count_gap_windows()))


  # load a sequence
  sequence = regions[0].get_region_as_rd_mean_sequences_with_windows(size=None,
                                                                   window_type=WindowType.from_string(configuration["HMM"]["train_windowtype"]),
                                                                   n_seqs=configuration["HMM"]["train_n_sequences_per_source"])

  observations = []
  for i in range(len(sequence)):
    observations.append(sequence[i][0])

  print("Sequence length: ",len(sequence))

  time_start = time.perf_counter()
  viterbi_path = hmm.viterbi(observations)
  time_end = time.perf_counter()
  print("Done. Execution time"
          " {0} secs".format(time_end - time_start))
  print("Log-probability of ML Viterbi path: ", viterbi_path[0])


  if viterbi_path[1] is not None:
    print("Viterbi path length: ", len(viterbi_path[1]))

    filename="viterbi_path.txt"
    counter = 0
    with open(filename, 'w') as f:
        f.write(str(len(viterbi_path[1])-1) + "\n")
        for item in range(len(sequence)):

            if sequence[item][0] == (-999.0, -999.0):
                counter += 1

            f.write(str(item)+ ":" + str(sequence[item][1]) + ":" + str(sequence[item][0]) + ":" + viterbi_path[1][item+1][1].name + "\n")
            #print("sequnce item: {0} state {1}".format(sequence[item], viterbi_path[1][item+1][1].name))
    print("There should be {0} gaps".format(counter))


if __name__ == '__main__':
  main()




