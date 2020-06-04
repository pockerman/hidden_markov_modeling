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
from cluster import Cluster
from region import Region

def compute_dbi(cluster_files, regionfile):

  # load the region
  region = Region.load(filename=regionfile)
  region.check_windows_sanity()
  region.get_mixed_windows()

  # load the clusters
  clusters=[]

  for file in cluster_files:
     cluster = Cluster.load(filename=file)
     cluster.windows = region.get_mixed_windows()
     clusters.append(cluster)


  dbi = Cluster.dbi(clusters)
  print("Dbi index is: {0}".format(dbi))


def main():


  path = "/home/a/ag568"


  print("Three clusters")
  clusters = ["cluster_0_MANHATAN_3_MEAN_RATIO.txt",
              "cluster_1_MANHATAN_3_MEAN_RATIO.txt",
              "cluster_2_MANHATAN_3_MEAN_RATIO.txt"]

  for cluster in clusters:
    cluster = path + "/" + cluster

  region = path + "/" + "region_0_MANHATAN_3_MEAN_RATIO.txt"
  compute_dbi(cluster_files=clusters, regionfile=region)


  print("Two clusters")
  clusters = ["cluster_0_MANHATAN_2_MEAN_RATIO.txt",
              "cluster_1_MANHATAN_2_MEAN_RATIO.txt",
              ]

  for cluster in clusters:
    cluster = path + "/" + cluster

  region = path + "/" + "region_0_MANHATAN_2_MEAN_RATIO.txt"
  compute_dbi(cluster_files=clusters, regionfile=region)

if __name__ == '__main__':

  main()

