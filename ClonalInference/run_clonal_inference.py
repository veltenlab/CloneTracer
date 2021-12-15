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
  parser.add_argument("--number_iterations", "-t", help="Number of iterations for inferring clonal hierarchies. Defaults to 300 for samples with only SNVs and 500 otherwise. Minimum 60.", 
                      type=int, default=300)
  parser.add_argument("--all_trees", "-a", help="Boolean indicating whether to run the model for all possible trees. Defaults to False which means a heuristic tree building approach is used.", 
                      type=bool, default=False)  

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
    
def create_tree_class(input_file, out_dir, name, single, cnv_celltype, gpu):
    
    # load data from JSON file
    with open(input_file, "rb") as f:
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
        if entry in input_data and input_data[entry]: 

            data_svi[entry] = torch.Tensor(input_data[entry])

        # otherwise set values to 0    
        else:
            data_svi[entry] = torch.zeros(nmuts)

    # priors for heteroplasmy (default values are 1000,1000 for nuclear, 1,1 for mitochondria and 2,100 for CNVs)
    for entry in ["h_alpha", "h_beta"]:

        # if present add information to dictionary
        if entry in input_data and input_data[entry]: 

            data_svi[entry] = torch.Tensor(input_data[entry])

        # otherwise set to default values
        else:

            h_mapper = {0: 2, 1: 1000, 2: 1}

            data_svi[entry] = torch.Tensor([h_mapper[mut] for mut in input_data["mut_type"]])  


    # add additional information for celltype-specific CNV model (if present)
    for entry in ["class_assign", "class_names", "celltype", "celltype_names", "cnv_ct_mean", "cnv_ct_sd", "umapx", "umapy"]:

        if entry in input_data and input_data[entry]: 

            if entry == "class_assign" and gpu:

                data_svi[entry] = torch.cuda.IntTensor(input_data[entry])

            elif entry == "class_assign" and not gpu:

                data_svi[entry] = torch.IntTensor(input_data[entry])

            elif entry == "cnv_ct_mean":

                data_svi[entry] = torch.Tensor(input_data[entry])

            elif entry in ["cnv_ct_sd", "celltype"]:

                data_svi[entry] = torch.tensor(input_data[entry])

            else:

                data_svi[entry] = input_data[entry]

        else:

            data_svi[entry] = []

    # rename bulk data entries
    data_svi["af_alpha"] = data_svi["bulk_M"]
    data_svi["af_beta"] = data_svi["bulk_N"]

    t = tree(name, data_svi)
    
    return(t)

# create tree class
t = create_tree_class(json_in, out_dir, name, single, cnv_celltype, gpu)

# set number of iterations (defaults to 300 when only SNVs are present and 500 when CNVs are present)
num_iter = args.number_iterations

if num_iter != 300:
    continue
else:
    if torch.nonzero(t.mut_type == 0).size()[0] == 0:
        
        num_iter = 300
        
    else:
        
        num_iter = 500
    
init = num_iter - 100

# run model
if not args.all_trees:

    # infer clonal hierarchy and compute clonal probabilities for the selected trees
    t.infer_hierarchy(num_iter, init, out_dir)
    
else:
    
    # run model for all possible trees
    t.run_all_trees(out_dir = out_dir, num_iter = num_iter)
