# executable python script to the clonal inference model in one sample

# Import packages
import os
import torch
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import pickle
import dill
import json
import argparse
import time
from datetime import datetime
from timeit import default_timer as timer
from datetime import timedelta

import pyro
import pyro.distributions as dist
import pyro.poutine as poutine
from pyro.infer import SVI, Trace_ELBO, TraceEnum_ELBO, config_enumerate, Predictive
from pyro.optim import Adam, AdagradRMSProp
from torch.distributions import constraints
pyro.set_rng_seed(100)

# parse command line specifications
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", "-i", help="Input JSON file. See a detailed explanation of the format in https://github.com/veltenlab/AMLScripts/tree/master/ClonalInference", type=str)
  parser.add_argument("--name", "-n", help="Sample name", type=str)
  parser.add_argument("--output_dir", "-o", help="Output directory. Results will be saved as a pickle file.", type=str)
  parser.add_argument("--single_sample", "-s", help="Boolean indicating whether the bulk data comes from one or more samples (e.g. different time points, T cells and myeloid cells)", 
                      type=bool, default=True)
  parser.add_argument("--cnv_celltype", "-c", help="Boolean indicating whether to use the celltype-specific model for CNVs", type=bool, default=False)
  parser.add_argument("--gpu", "-g", help="Boolean indicating whether to use a GPU for model inference", type=bool, default=False)  
  args = parser.parse_args()

# import helper functions to run the model
from helper_functions import *

# load parse arguments
json_in = args.input
out_dir = args.output_dir
name = args.name
single = args.single_sample
cnv_celltype = args.cnv_celltype
gpu = args.gpu

# set DoubleTensor as default and specify whether to use GPU
if gpu:
    torch.set_default_tensor_type(torch.cuda.DoubleTensor)
else:
    torch.set_default_tensor_type(torch.DoubleTensor)
    
# load data from JSON file
with open(json_in, "rb") as f:
    input_data = json.load(f)
    
# get number of mutations
nmuts = len(input_data["mut_names"])

# if there are more than one sample 
    
# make a dictionary with input data    
data_svi = {"M": torch.Tensor(input_data["M"]),
             "N": torch.Tensor(input_data["N"]),
             "mut_type": torch.Tensor(input_data["mut_type"]),
             "names": input_data["mut_names"],
             "barcodes": input_data["cell_barcode"],
             "class_af": single,
             "cnv_celltype": cnv_celltype}

# bulk data if present
for entry in ["bulk_M", "bulk_N", "r_cnv"]:
    
    # if present add information to dictionary
    if entry in input_data and input_data["entry"]: 
    
        data_svi[entry] = torch.Tensor(input_data[entry])
        
    # otherwise set values to 0    
    else:
        data_svi[entry] = torch.zeros(nmuts)
        
# priors for heteroplasmy (default values are 1000,1000 for nuclear, 1,1 for mitochondria and 2,100 for CNVs)
for entry in ["h_alpha", "h_beta"]:
    
    # if present add information to dictionary
    if entry in input_data and input_data["entry"]: 
    
        data_svi[entry] = torch.Tensor(input_data[entry])
        
    # otherwise set to default values
    else:
        
        h_mapper = {0: 2, 1: 1000, 2: 1}
        
        data_svi[entry] = torch.Tensor([h_mapper[mut] for mut in input_data["mut_type"]])  
        
        
# add additional information for celltype-specific CNV model (if present)
for entry in ["class_assign", "class_names", "celltype", "celltype_names", "cnv_ct_mean", "cnv_ct_sd", "umapx", "umapy"]:
    
    if entry in input_data and input_data["entry"]: 
        
        if entry in ["class_assign", "celltype"] and gpu:
            
            data_svi[entry] = torch.cuda.IntTensor(input_data[entry])
        
        elif entry in ["class_assign", "celltype"] and not gpu:
            
            data_svi[entry] = torch.IntTensor(input_data[entry])
            
        elif entry == "cnv_ct_mean":
            
            data_svi[entry] = torch.Tensor(input_data[entry])
            
        elif entry == "cnv_ct_sd":
            
            data_svi[entry] = torch.tensor(input_data[entry])
            
        else:
            
            data_svi[entry] = input_data[entry]
            
    else:
        
        data_svi[entry] = []
            

with open(out_dir+"/test_dict.pickle") as t:
    pickle.dump(data_svi, t, protocol=pickle.HIGHEST_PROTOCOL)
