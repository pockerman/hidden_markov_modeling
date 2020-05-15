# Hidden Markov Model

## Dependencies

- ```pyclustering```
- ```pysam```
- ```pomegranate```

## How to use

Execute the ```train.py``` script:

```
python3 train.py
```

The ```train.py``` reads the local ```config.json``` file
that provides configuration parameters for the training.
Below is an explanation of the entries in ```config.json```

- ```regions```: is a map with ```start``` and ```end```
lists of the regions to be used. Regions do not have to be consecutive
```
 "regions":{
      "start":[1000000],
      "end":[10000000]
  },
```

- ```"chromosome"```: is single valued and specifies which chromosome to use

```
"chromosome": "chr1",
```

- ```reference_file```: specifies the refernce sequence file
- ```no_wga_file```: specifies the non-WGA-treated sequence file
- ```test_file```: specifies the WGA-treated sequence file
- ```debug```: a map with debug options currently only ```log_bam_for_debug```
flag is supported. If set to ```true``` the bam entries read for the regions are written
in the ```logger_file``` .
- ```"window_size"```: The window size to use
- ```"fill_missing_window_data"```: flag indcating if missing window RD counts should be filled
- ```"fill_missing_window_data_factor"```: The factor to use for filling missing window data
- ```"quality_threshold"```: Threshold to read a sequence from bam. If ```null`` no threshold is used,
- ```"save_windows"```: Flag indicating if the created windows should be saved,
- ```"windows_filename"```:The fille to save the windows,
- ```"remove_windows_with_N"```: Flag indicating that any window containing an ```N``` should be removed,
- ```"mark_N_windows"```: Flag indicating whethre to mark windows containig an ```N```,
- ```"mark_for_N_windows"```: The mark to use for the ```N``` windows
-  ```"outlier_remove"```: Map indicating how to mark outlier windows
```
 "outlier_remove":{

    "name":"zscore",
    "config":{
      "sigma_factor":2
    }
  },
```
- ```"clusterer"```: A map with the properties of the clustering method to be used

```
 "clusterer":{
    "name":"kmedoids",
    "config":{
    "init_cluster_idx":"random_from_data", or a list of indices
    "metric":"MANHATAN", or "EUCLIDEAN"
    "features":["mean"]
    }
  },
```

- ```"label_clusters"```: Flag indicating if the clusters should be labelled
- ```"labeler"```: If ```label_clusters``` is set to ```true``` this map contains
the properties of the labeler

```
"label_clusters": true,
  "labeler":{
    "name":"mean_diff",
    "tuf_mean_min": 1.5,
    "tuf_mean_max": 8.5,
    "states":["DELETE", "OTHER", "OTHER", "TUF"]
  },
```

- ```"cluster_distribution"```: Map indicating how the densities of each formed cluster should be modelled

```
"cluster_distribution":{
    "name": "gmm",
    "config":{
      "distributions":{
          "TUF":{
              "type":"gmm",
              "dists":["normal", "uniform"],
              "weights":null,
              "uniform_params":[2.0, 8.0] # U(2.0, 8.0)
          },
          "OTHER":{
              "type":"distribution", # use gmm for multinomial
              "dists":["normal"]
          },
          "DELETE":{
              "type":"distribution",
              "dists":["normal"]
          }
      }
    }
  },
```
-   ```"save_cluster_densities"```: Flag indicating if the cluster data should be saved,
- ```"HMM"```: Properties for HMM model

```
"HMM": {
    "name":"HMM_Model", # name of the model
    "train":false, # whether train or not
    "train_solver": "baum-welch", # if train how: "baum-welch" or "viterbi"
    "lr_decay":0.7,
    "inertia":null,
    "verbose":true,
    "save_model":true,
    "save_hmm_filename":"HMM_Model",
    "start_prob":{
      "OTHER_0": 0.475, # not used
      "OTHER_1": 0.475, # not used
      "DELETE": 0.0166,
      "TUF": 0.0166
    },
    "train_sequence_size": 100,
    "train_sequence_source": "region",
    "train_n_sequences_per_source":10000,
    "train_windowtype":"both"
  },
```

- ```"logger_file"```: The filename for the logs,
- ```"logger_level"```: Minimum level of logging








