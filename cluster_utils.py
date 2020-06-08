import sys
import numpy as np
from pomegranate import *


from exceptions import Error
from helpers import INFO, WARNING
from helpers import timefn
from helpers import WindowType,  WindowState
from preprocess_utils import get_distributions_list_from_names


def similar_to_cluster(data, clusters):

  cluster_ids = []

  for item in data:

    dist = sys.float_info.max
    cidx = -1
    for cluster in clusters:

      dist_clst = cluster.distance_from(item)

      if dist_clst < dist:
        dist = dist_clst
        cidx = cluster.cidx

    cluster_ids.append((item, cidx))
  return cluster_ids


def build_cluster_mean_and_std(clusters, **kwargs):

    for cluster in clusters:

      indeces = cluster.indexes

      wga_data = np.empty((1,0), float)
      no_wga_data = np.empty((1,0), float)
      windows = cluster.windows

      for idx in indeces:
        window = windows[idx]
        if window.is_gap_window() == False:
          mu1, mu2 = window.get_rd_statistic(statistics="mean",
                                             name=WindowType.BOTH)

          wga_data = np.append(wga_data, np.array(mu1))
          no_wga_data = np.append(no_wga_data, np.array(mu2))


      # collected the data create the GMM for each
      # component in the cluster
      wga_params={"mean": np.mean(wga_data),
                  "std": np.std(wga_data)}

      no_wga_params={"mean": np.mean(no_wga_data),
                     "std": np.std(no_wga_data)}

      cluster.wga_mean = np.mean(wga_data)
      cluster.wga_std = np.std(wga_data)
      cluster.no_wga_mean = np.mean(no_wga_data)
      cluster.no_wga_std = np.std(no_wga_data)

def find_name_from_state(state, **kwargs):

  for clst in kwargs['clusters']:
    if kwargs['clusters'][clst]['state'] == state:
      return clst

  raise Error("No cluster with state {0} found".format(state))


@timefn
def build_cluster_densities(clusters_lst, **kwargs):

      print("{0} Build cluster densities".format(INFO))

      for cluster in clusters_lst:

        if isinstance(cluster.state, str):
          state = cluster.state.lower()
        else:
          state = cluster.state.name.lower()

        print("{0} Cluster state: {1}".format(INFO, state))
        name = find_name_from_state(state, **kwargs)
        print("{0} Cluster name: {1}".format(INFO, name))

        # collected the data create the GMM for each
        # component in the cluster
        wga_params={"mean": cluster.wga_mean,
                    "std": cluster.wga_std}

        no_wga_params={"mean": cluster.no_wga_mean,
                       "std": cluster.no_wga_std}

        if state == 'tuf':

          if 'names' in kwargs['clusters'][name]["distributions"]["wga"] and \
            'uniform' in kwargs['clusters'][name]["distributions"]["wga"]['names']:
            uniform_params = kwargs['clusters'][name]["distributions"]["wga"]["uniform"]["params"]
            wga_params["uniform_params"] = uniform_params

          if 'names' in kwargs['clusters'][name]["distributions"]["no_wga"] and \
            'uniform' in kwargs['clusters'][name]["distributions"]["no_wga"]['names']:
              uniform_params = kwargs['clusters'][name]["distributions"]["no_wga"]["uniform"]["params"]
              no_wga_params["uniform_params"] = uniform_params
        else:
          if 'names' in kwargs['clusters'][name]["distributions"]["wga"] and \
          'uniform' in kwargs['clusters'][name]["distributions"]["wga"]['names']:
            uniform_params = kwargs['clusters'][name]["distributions"]["wga"]["uniform"]["params"]
            wga_params["uniform_params"] = uniform_params
          elif 'name' in kwargs['clusters'][name]["distributions"]["wga"] and \
            kwargs['clusters'][name]["distributions"]["wga"]['name'] == 'uniform':
            uniform_params = kwargs['clusters'][name]["distributions"]["wga"]["uniform"]["params"]
            wga_params["uniform_params"] = uniform_params

          if 'names' in kwargs['clusters'][name]["distributions"]["no_wga"] and \
            'uniform' in kwargs['clusters'][name]["distributions"]["no_wga"]['names']:
              uniform_params = kwargs['clusters'][name]["distributions"]["no_wga"]["uniform"]["params"]
              no_wga_params["uniform_params"] = uniform_params
          elif 'name' in kwargs['clusters'][name]["distributions"]["no_wga"] and \
            kwargs['clusters'][name]["distributions"]["no_wga"]['name'] == 'uniform':
              uniform_params = kwargs['clusters'][name]["distributions"]["no_wga"]["uniform"]["params"]
              no_wga_params["uniform_params"] = uniform_params


        type_ = kwargs['clusters'][name]["distributions"]["wga"]["type"]
        if type_ == "gmm":

          names = kwargs['clusters'][name]["distributions"]["wga"]["names"]
          weights = kwargs['clusters'][name]["distributions"]["wga"]["weights"]

          wga_gmm = \
              GeneralMixtureModel(
                get_distributions_list_from_names(names,
                                                  wga_params),
                                  weights=weights)

          cluster.wga_density = wga_gmm
        elif type_ == 'distribution':
          dist_name = kwargs['clusters'][name]["distributions"]["wga"]["name"]


          if dist_name == 'uniform':
            wga_params['uniform_params'] = \
              kwargs['clusters'][name]["distributions"]["wga"]["parameters"]

          wga_dist = get_distributions_list_from_names([dist_name],
                                                       wga_params)[0]
          cluster.wga_density = wga_dist
        elif type_ is None:
          print("{0} No density specified for WGA sample".format(WARNING))
        else:
          raise Error("Invalid cluster "
                      "distribution method for WGA sample."
                      " Method: {0}".format(type_))


        type_ = kwargs['clusters'][name]["distributions"]["no_wga"]["type"]
        if type_ == "gmm":

          names = kwargs['clusters'][name]["distributions"]["no_wga"]["names"]
          weights = kwargs['clusters'][name]["distributions"]["no_wga"]["weights"]
          non_wga_density = \
              GeneralMixtureModel(get_distributions_list_from_names(names,
                                                                    no_wga_params),
                                  weights=weights)

          cluster.no_wga_density = non_wga_density

        elif type_ == "distribution":
          dist_name = kwargs['clusters'][name]["distributions"]["no_wga"]["name"]

          if dist_name == 'uniform':
            no_wga_params['uniform_params'] = \
              kwargs['clusters'][name]["distributions"]["wga"]["parameters"]

          non_wga_density = get_distributions_list_from_names([dist_name],
                                                              no_wga_params)[0]
          cluster.no_wga_density = non_wga_density
        elif type_ is None:
          print("{0} No density specified for NO_WGA sample".format(WARNING))

        else:
            raise Error("Invalid cluster "
                        "distribution method for NO_WGA sample. "
                        " Method: {0}".format(type_))

      return clusters_lst

