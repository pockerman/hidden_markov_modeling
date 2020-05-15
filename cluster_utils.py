import numpy as np
from sklearn.neighbors import KernelDensity
from pomegranate import *
from exceptions import Error
from helpers import INFO
from helpers import WindowType,  WindowState


def label_clusters(clusters, method, **kwargs):
    
    labeler = kwargs["labeler"]
    if method == "mean_diff":
        
        print("{0} Labeler: {1}".format(INFO, method))
        
        for cluster in clusters:
            
            print("{0} Cluster {1}".format(INFO, cluster.cidx))
            
            mu_wga = cluster.get_statistics(statistic="mean",
                                            window_type=WindowType.WGA)
            
            
            mu_no_wga = cluster.get_statistics(statistic="mean",
                                               window_type=WindowType.NO_WGA)
            
            print("{0} Cluster WGA mean: {1}".format(INFO, mu_wga))
            print("{0} Cluster NO WGA mean: {1}".format(INFO, mu_no_wga))
            
            if mu_wga >= labeler["tuf_mean_min"] and mu_wga <= labeler["tuf_mean_max"]:
                
                # this is potential tuf
                if np.fabs(mu_no_wga - mu_wga) > 5.0:
                    cluster.state = WindowState.TUF
                else:
                    cluster.state = WindowState.OTHER
                    
            elif mu_wga < labeler["tuf_mean_min"] and mu_no_wga < labeler["tuf_mean_min"]:
                cluster.state = WindowState.DELETE
                
            else:
                cluster.state = WindowState.OTHER
                
            print("{0} Cluster state is {1}".format(INFO, cluster.state.name))

    return clusters        

    

def build_cluster_densities(clusters, **kwargs):
  """
  Establish the probability destributions underlying
  the data for each cluster

  Parameters
  ----------
  clusters : TYPE
    DESCRIPTION.
  **kwargs : TYPE
    DESCRIPTION.

  Returns
  -------
  None.

  """
  print("{0} Type of cluster density fitted: {1}".format(INFO,
                                                         kwargs["name"]))
  
  config = kwargs["config"]
  
  if kwargs["name"] == "gmm":

    # distributions for the first component
    #wga_dist = kwargs["config"]["wga_dist"]
    #wga_weights = kwargs["config"]["wga_weights"]

    #print("{0} Distributions for WGA {1} ".format(INFO, wga_dist))

    #non_wga_dist = kwargs["config"]["no_wga_dist"]
    #no_wga_weights = kwargs["config"]["no_wga_dist_weights"]

    #print("{0} Distributions for NON WGA {1} ".format(INFO, wga_dist))

    for cluster in clusters:
        
      clust_dists = config["distributions"][cluster.state.name]
      indeces = cluster.indexes

      wga_data = np.empty((1,0), float)
      no_wga_data = np.empty((1,0), float)
      windows = cluster.windows

      for idx in indeces:
        window = windows[idx]
        if window.is_n_window() == False:
          mu1, mu2 = window.get_rd_stats(statistics="mean",
                                         name=WindowType.BOTH)

          wga_data = np.append(wga_data, np.array(mu1))
          no_wga_data = np.append(no_wga_data, np.array(mu2))
          
          
      # collected the data create the GMM for each
      # component in the cluster
      wga_params={"mean": np.mean(wga_data),
                  "std": np.std(wga_data)}
      
      no_wga_params={"mean": np.mean(no_wga_data),
                     "std": np.std(no_wga_data)}
      
      if cluster.state.name == 'TUF':
          wga_params["uniform_params"] = config["distributions"]["TUF"]["uniform_params"]
          no_wga_params["uniform_params"] = config["distributions"]["TUF"]["uniform_params"]
      
      if config["distributions"][cluster.state.name]["type"] == "gmm":
          wga_gmm = \
              GeneralMixtureModel(_get_distributions_list_from_names(clust_dists["dists"],
                                                                     wga_params),
                                  weights=config["distributions"][cluster.state.name]["weights"])

          cluster.wga_density = wga_gmm

          non_wga_density = \
              GeneralMixtureModel(_get_distributions_list_from_names(clust_dists["dists"],
                                                               no_wga_params),
                                  weights=config["distributions"][cluster.state.name]["weights"] )

          cluster.no_wga_density = non_wga_density
          
      elif config["distributions"][cluster.state.name]["type"] == "distribution":
         wga_dist = _get_distributions_list_from_names(clust_dists["dists"], 
                                                       wga_params)[0]
         cluster.wga_density = wga_dist
         
         non_wga_density = _get_distributions_list_from_names(clust_dists["dists"], 
                                                              no_wga_params)[0]
         cluster.no_wga_density = non_wga_density
                                  
  else:
    raise Error("Invalid cluster distribution method")

  return clusters


def save_clusters_desnity(clusters, **kwargs):

  for cluster in clusters:
    filename = "cluster_" + str(cluster.cidx) + "_density.txt"
    save_cluster_density(cluster=cluster, filename=filename, **kwargs)


def clusters_statistics(clusters):
  """
  Claculate various statistics for the windows
  clustered in clusters

  Parameters
  ----------
  clusters : list of lists
    Contains the clustered window indexes. There are
    len(clusters) clusters
  windows : list of Window objects
    DESCRIPTION.

  Returns
  -------
  a map with the calculated statistics for each
  cluster

  """

  statistics = {}

  for c in range(len(clusters)):

    statistics[c] = clusters[c].get_statistics(statistic="all")
  return statistics


def save_cluster_density(cluster, filename, **kwargs):


    if kwargs["name"] == "kde":

      arr = _form_cluster_2d_array(cluster=cluster)
      log_probs = cluster.density.score_samples(arr)

      with open(filename, 'w') as file:
        file.write(str(log_probs))
      return

    raise Error("Invalid cluster distribution method")


def _form_cluster_2d_array(cluster):
  indeces = cluster.indexes
  arr = np.empty((0, 2), float)
  windows = cluster.windows
  for idx in indeces:
    window = windows[idx]
    mu1, mu2 = window.get_rd_stats(statistics="mean", name="both")
    arr = np.append(arr, np.array([[mu1, mu2]]), axis=0)
  return arr


def _get_distributions_list_from_names(dists_name, params):

  dists = []

  for name in dists_name:
    if name == "normal":
      dists.append(NormalDistribution(params["mean"], params["std"]))
    elif name == "poisson":
      dists.append(PoissonDistribution(params["mean"]))
    elif name == "uniform":
      dists.append(UniformDistribution(params["uniform_params"][0] ,
                                       params["uniform_params"][1]))
    else:
      raise Error("Name {0} is an unknown distribution ".format(name))
  return dists

