# script to run the model on all possible trees of patient K

# Import packages
import os
import torch
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
import itertools
import os
import pickle

import dill

import time
from datetime import datetime
from timeit import default_timer as timer
from datetime import timedelta

import scipy
from scipy.stats import betabinom

import pyro
import pyro.distributions as dist
import pyro.poutine as poutine
from pyro.infer.autoguide import AutoDelta
from pyro.infer import SVI, Trace_ELBO, MCMC, NUTS, TraceEnum_ELBO, config_enumerate, infer_discrete
from pyro.optim import Adam, AdagradRMSProp, MultiStepLR
from torch.distributions import constraints
from pyro.ops.indexing import Vindex
pyro.set_rng_seed(100)

# import helper functions and tree class form source code files (present in the same directory)
from helper_functions import *

torch.set_default_tensor_type(torch.cuda.DoubleTensor)

# load data from patient K
with open("data/AKLK.json") as f:
    data_k = json.load(f)
    
# I will remove RUNX1 and the trisomy in chr8 for the moment
data_k["M"] = np.array(data_k["M"])[:, 1:6]
data_k["N"] = np.array(data_k["N"])[:, 1:6]
data_k["colnames"] = data_k["colnames"][1:6]
    
# add type of mutations: 0 = CNV, 1 = nuclear, 2 = mitochondrial
data_k["type"] = np.array([1, 1, 2, 2, 2])

# create tree class
data = {"M": torch.Tensor(data_k["M"]),
         "N": torch.Tensor(data_k["N"]),
         "mut_type": torch.Tensor([1,1,2,2,2]),
         "h_alpha": torch.Tensor([1000.0, 1000.0, 1, 1, 1]),
         "h_beta": torch.Tensor([1000.0, 1000.0, 1, 1, 1]),
         "af_alpha": torch.Tensor([[44, 43, 1, 1, 1], [0,0,0,0,0]]),
         "af_beta": torch.Tensor([[100-44, 100-43, 1, 1,1], [1,1,1,1,1]]),
         "r_cnv": torch.Tensor([0,0,0,0,0]),
         "names": data_k["colnames"], 
         "barcodes": data_k["cell_barcode"], 
         "umapx": data_k["umapx"], 
         "umapy": data_k["umapy"],
         "class_af": True, 
         "class_assign": torch.cuda.IntTensor(data_k["timepoint"]),
         "class_names": ["day0", "day15"], 
         "cnv_celltype": False,                                         
         "celltype": [],
         "celltype_names": [],
         "cnv_ct_mean": [],
         "cnv_ct_sd": []}

k = tree("K", data)

k.run_all_trees(outfile = "svi_objects/all_trees_K.pickle")