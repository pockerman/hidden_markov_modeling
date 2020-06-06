import time
import sys
import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from helpers import read_configuration_file
from helpers import set_up_logger
from helpers import timefn
from helpers import INFO
from cluster import Cluster
from region import Region


def load_region_mean(filename):

  with open(filename) as file:
        context = file.read()
        size = len(context)
        arraystr= context[1:size-1]
        arraystr = arraystr.split(',')
        region_means = [float(item) for item in arraystr]
        return region_means


def plot_hist(data, nbins=35, kde=False, rug=True):

  sns.distplot(data, bins=nbins, kde=kde, rug=rug)
  plt.show()

def plot_2d_hist(data1, data2, nbins=35, density=False, min_=0.0, max_=50.0):

  plt.hist2d(data1, data2,
             bins=[nbins, nbins], cmap='Blues', density=density,
             cmax=1000, cmin=0, alpha=0.99,
             range=((min_, max_), (min_, max_)))
  #plt.title = "Hist 2D"
  plt.show()


def main(config):

  wga_region_file = config["wga_region_file"]
  no_wga_region_file = config["no_wga_region_file"]

  wga_region_mean = load_region_mean(filename=wga_region_file)

  print("{0} WGA mean: {1}".format(INFO, np.mean(wga_region_mean)))
  print("{0} WGA var: {1}".format(INFO, np.var(wga_region_mean)))

  no_wga_region_mean = load_region_mean(filename=no_wga_region_file)

  print("{0} NO-WGA mean: {1}".format(INFO, np.mean(no_wga_region_mean)))
  print("{0} NO-WGA var: {1}".format(INFO, np.var(no_wga_region_mean)))

  X = np.stack((no_wga_region_mean, wga_region_mean), axis=0)
  print("{0} Covariance WGA-NO-WGA: {1}".format(INFO, np.cov(X)))


  # how much the WGA explains the NO-WGA

  plot_2d_hist(data1=no_wga_region_mean,
               data2=wga_region_mean,nbins=80)








if __name__ == '__main__':
    print("{0} Start clustering...".format(INFO))
    total_start = time.perf_counter()
    description = "Check the README file for "
    "information on how to use the script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=str, default='eda_config.json',
                        help="You must specify a json "
                        "formatted configuration file")


    print("{0} Read configuration file".format(INFO))
    args = parser.parse_args()
    configuration = read_configuration_file(args.config)
    print("{0} Done...".format(INFO))
    sys.stdout.flush()

    main(config=configuration)
    total_end = time.perf_counter()
    print("{0} Finished clustering. "
          "Total execution time {1} secs".format(INFO, total_end - total_start))
    sys.stdout.flush()
