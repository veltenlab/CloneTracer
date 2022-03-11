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
  parser.add_argument("--input", "-i", help="Input JSON file. See a detailed explanation of the format in https://github.com/veltenlab/AMLScripts/tree/master/ClonalInference", 
                      type=str)
  parser.add_argument("--name", "-n", help="Sample name", type=str)
  parser.add_argument("--output_dir", "-o", help="Output directory. Results will be saved as a pickle file.", type=str)
  parser.add_argument("--multiple_samples", "-s", 
                      help="Boolean indicating whether the bulk data comes from more than one sample (e.g. different time points, T cells and myeloid cells). Default: False", 
                      default=False, action="store_true")
  parser.add_argument("--cnv_celltype", "-c", help="Boolean indicating whether to use the celltype-specific model for CNVs", default=False, action="store_true")
  parser.add_argument("--gpu", "-g", help="Boolean indicating whether to use a GPU for model inference", default=False, action="store_true")   
  parser.add_argument("--number_iterations", "-t", 
                      help="Number of iterations for inferring clonal hierarchies. Defaults to 300 for samples with only SNVs and 500 otherwise. Minimum 60.", 
                      type=int, default=0)
  parser.add_argument("--all_trees", "-a", 
                      help="Boolean indicating whether to run the model for all possible trees. Defaults to False which means a heuristic tree building approach is used.", 
                      default=False, action="store_true")  

args = parser.parse_args()

# import helper functions to run the model
from helper_functions import *

# load parse arguments
json_in = args.input
out_dir = args.output_dir
name = args.name
mult_samp = args.multiple_samples
cnv_celltype = args.cnv_celltype
gpu = args.gpu

# create output directory if it does not exist
if not os.path.isdir(out_dir):
    
    os.mkdir(out_dir)

# set DoubleTensor as default and specify whether to use GPU
if gpu:
    torch.set_default_tensor_type(torch.cuda.DoubleTensor)
else:
    torch.set_default_tensor_type(torch.DoubleTensor)
    
# create tree class
t = create_tree_class(json_in, name, mult_samp, cnv_celltype, gpu)

# set number of iterations (defaults to 300 when only SNVs are present and 500 when CNVs are present)
num_iter = args.number_iterations

if num_iter == 0:
    
    # only SNVs
    if torch.nonzero(t.mut_type == 0).size()[0] == 0:
        
        num_iter = 300
        
    # sample with CNVs
    else:
        
        num_iter = 500
    
init = num_iter - 100

# run model
if args.all_trees == False:

    # infer clonal hierarchy and compute clonal probabilities for the selected trees
    t.infer_hierarchy(num_iter, init, out_dir)
    
else:
    
    # run model for all possible trees
    t.run_all_trees(out_dir = out_dir, num_iter = num_iter)
