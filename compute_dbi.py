import time
import sys

from helpers import timefn
from helpers import INFO
from cluster import Cluster
from region import Region

@timefn
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


  for cluster in clusters:
     cluster.diameter

     for other in clusters:
       cluster.calculate_distance_from_other(other=other)

  for cluster in clusters:
     print("{0} Cluster {1}".format(INFO, cluster.cidx))
     print("{0} Cluster diameter {1}".format(INFO, cluster.diameter))
     print("{0} Cluster other distance {1}".format(INFO, cluster.distance_from_others))


  dbi = Cluster.dbi(clusters)
  print("Dbi index is: {0}".format(dbi))


def main():


  path = "/home/a/ag568"


  print("Three clusters")
  clusters = ["cluster_0_MANHATAN_3_MEAN_RATIO.txt",
              "cluster_1_MANHATAN_3_MEAN_RATIO.txt",
              "cluster_2_MANHATAN_3_MEAN_RATIO.txt"]

  for clst in range(len(clusters)):
    clusters[clst] = path + "/" + clusters[clst]

  region = path + "/" + "region_0_MANHATAN_3_MEAN_RATIO.txt"
  compute_dbi(cluster_files=clusters, regionfile=region)


  print("Two clusters")
  clusters = ["cluster_0_MANHATAN_2_MEAN_RATIO.txt",
              "cluster_1_MANHATAN_2_MEAN_RATIO.txt",
              ]

  for clst in range(len(clusters)):
    clusters[clst] = path + "/" + clusters[clst]

  region = path + "/" + "region_0_MANHATAN_2_MEAN_RATIO.txt"
  compute_dbi(cluster_files=clusters, regionfile=region)

if __name__ == '__main__':

  main()

