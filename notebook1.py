import json
from pomegranate import*
from hmm_helpers import build_hmm



def main():
  sample = "sample_1"
  ffile = "../"
  hmm_file = ffile + "HMM_Model_0.json"


  #import pdb
  #pdb.set_trace()
  hmm = build_hmm(hmm_file)
  json_str = hmm.to_json()

  print(json_str)

  #new_hmm = HiddenMarkovModel.from_json(json_str)


if __name__ == '__main__':
  main()
