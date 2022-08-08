# file with helper functions to find the best tree(s) heuristically

# Import packages
import os
import sys
import pandas as pd
import torch
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json
import math
import itertools
import os
import pickle
import copy

import dill

import time
from datetime import datetime
from timeit import default_timer as timer
from datetime import timedelta


import pyro
import pyro.distributions as dist
import pyro.poutine as poutine
from pyro.infer import SVI, Trace_ELBO, TraceEnum_ELBO, config_enumerate
from pyro.optim import Adam, AdagradRMSProp, MultiStepLR
from torch.distributions import constraints
pyro.set_rng_seed(100)

# function to create tree class object from JSON input file
def create_tree_class(input_file, name, mult_samp, cnv_celltype, gpu):
    

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
                 "class_af": mult_samp,
                 "cnv_celltype": cnv_celltype}

    # bulk data if present
    for entry in ["bulk_M", "bulk_N", "r_cnv"]:

        # if present add information to dictionary
        if entry in input_data and input_data[entry]:

            data_svi[entry] = torch.Tensor(input_data[entry])
            
            if entry == "bulk_M":
                
                data_svi["bulk_af"] = True

        # otherwise set values to 0
        else:
            data_svi[entry] = torch.zeros(nmuts)
            
            if entry == "bulk_M":
                
                data_svi["bulk_af"] = False

    # priors for heteroplasmy (default values are 1000,1000 for nuclear, 1,1 for mitochondria and 2,100 for CNVs)
    for entry in ["h_alpha", "h_beta"]:

        # if present add information to dictionary
        if entry in input_data and input_data[entry]:

            data_svi[entry] = torch.Tensor(input_data[entry])

        # otherwise set to default values
        else:

            if entry == "h_alpha":

                h_mapper = {0: 2, 1: 1000, 2: 1}

                data_svi[entry] = torch.Tensor([h_mapper[mut] for mut in input_data["mut_type"]])

            else:

                h_mapper = {0: 100, 1: 1000, 2: 1}

                data_svi[entry] = torch.Tensor([h_mapper[mut] for mut in input_data["mut_type"]])
        
    # convert dictionary with cnv priors to tensor    
    if cnv_celltype:
        
        mean = list()
        sd = list()
        
        for chrom in input_data["cnv_priors"]:
            
            mean.append(input_data["cnv_priors"][chrom]["mean"])
            sd.append(input_data["cnv_priors"][chrom]["sd"])
            
        data_svi["cnv_ct_mean"] = torch.Tensor(mean)
        data_svi["cnv_ct_sd"] = torch.Tensor(sd)
        
    else:
        
        data_svi["cnv_ct_mean"] = []
        data_svi["cnv_ct_sd"] = []

    # add additional information for celltype-specific CNV model (if present)
    for entry in ["class_assign", "class_names", "celltype", "celltype_names", "umapx", "umapy"]:

        if entry in input_data and input_data[entry]:

            if entry in ["class_assign"] and gpu:

                data_svi[entry] = torch.cuda.IntTensor(input_data[entry])

            elif entry in ["class_assign"] and not gpu:

                data_svi[entry] = torch.IntTensor(input_data[entry])
                
            elif entry in ["celltype"]:
                
                data_svi[entry] = torch.tensor(input_data[entry])

#             elif entry == "cnv_ct_mean":

#                 data_svi[entry] = torch.Tensor(input_data[entry])

#             elif entry in ["cnv_ct_sd", "celltype"]:

#                 data_svi[entry] = torch.tensor(input_data[entry])

            else:

                data_svi[entry] = input_data[entry]

        else:

            data_svi[entry] = []

    # rename bulk data entries
    data_svi["af_alpha"] = data_svi["bulk_M"]
    data_svi["af_beta"] = data_svi["bulk_N"]

    t = tree(name, data_svi)

    return(t)


# function to get parameters after svi
def get_params(svi):
    svi_params = [{i: k.cpu().detach().numpy()} for i,k in pyro.get_param_store().items()]
    params = {"svi_params": {}}
    for i in range(len(svi_params)):
        params["svi_params"].update(svi_params[i])
    return(params)

# get allele frequencies from clonal fractions
def get_af(cf, A):
    noroot = torch.transpose(A,0,1)
    return(torch.matmul(noroot, cf)/2)

# function to make a boolean torch based on node and mutation
def bool_node(tree, j, k):
    
    list_bool = []
    for node in k:
        if tree[node,j] == 1:
            list_bool.append(True)
            
        else:
            list_bool.append(False)
            
    return(torch.tensor(list_bool))


# define a class tree 
class tree:
    # instance attribute
    def __init__(self, name, data):
        self.name = name
        self.tree = []
        # add potential trees 
        self.potential_trees = [torch.Tensor([[0,0],
                                              [1,0],
                                              [1,1]]),
                                torch.Tensor([[0,0],
                                              [0,1],
                                              [1,1]]),
                                torch.Tensor([[0,0],
                                              [1,0],
                                              [0,1]]),
                                torch.Tensor([[0,0],
                                              [1,1]])]
        # mutations included in the tree
        self.muts = []
        # list of children for each tree
        self.children = [[[1,2],[2],[]],[[1,2],[2],[]], [[1,2],[],[]], [[1],[]]]
        # list of parents for each tree
        self.parents = [[[], [0], [0,1]], [[], [0], [0,1]], [[], [0], [0]], [[], [0]]]
        # attribute to store the ELBO during SVI
        self.elbo = []
        # attribute to store the iteration with lowest ELBO for each tree
        self.best_iter = []
        # attribute to store the parameter values of the tree guides
        self.parameters = []
        # attribute to store the indices of selected trees
        self.tree_indices = []
        # store number of iterations for SVI
        self.num_iter = 0
        # store generated posterior distributions
        self.posterior = []
        # store final parameters of posterior distribution
        self.params_post = []
        # store clonal assignment posterior for each selected tree
        self.clone_probs = []
        # store posterior predictive distributions of latent variables
        self.post_predictive = {}
        # store posterior predictive distributions for post. predictive checks 
        self.post_checks = {}
        
        for key, value in data.items():
            setattr(self, key, value)
    
    # select mutations with highest coverage
    def sel_mutation(self, initial = False):
        
        # get mutations not included in the tree
        nmut = list(set([*range(self.M.shape[1])])-set(self.muts))
        
        if len(nmut) == 0:
            print("All mutations have been added.")
            return(0)
        
        # if there is only one mutation left include it
        if len(nmut) == 1:
            self.muts.append(nmut[0])
            return(nmut[0])
        else:
            zeros = torch.zeros(len(nmut))
            # for CNVs only the counts in M matrix are considered
            t = torch.where(self.mut_type > 1., 1., self.mut_type)
            # iterate through mutations
            for ind,c in enumerate(nmut):

                # sum the total counts
                counts = self.M + self.N*t

                # count number of dropout cells 
                zeros[ind] = counts[:,c][counts[:,c] == 0].shape[0]

            if initial == True:
                # select top 2 mutations by covered cells     
                selected = torch.argsort(zeros)[0:2]
                for i in selected.tolist():
                    self.muts.append(i)

            else: 
                selected = nmut[torch.argsort(zeros)[0]]
                self.muts.append(selected)

            return(selected)
    
    #define model for SVI
    @config_enumerate
    def model(self, M, N, A, mut_type, h_alpha, h_beta, af_alpha, af_beta, names, r_cnv, 
              bulk_af = True, class_af = False, class_names = [], class_assign = [], 
              cnv_celltype = False, celltype = [], celltype_names = [], cnv_ct_mean = [], cnv_ct_sd = []):

        types = ["CNV", "nuclSNV", "mitoSNV"]

        nmut = M.shape[1]
        ncells = M.shape[0]
        
        # get number of mutations of each type
        ncnv = int(torch.sum(self.mut_type == 0.))
        nnucl = int(torch.sum(self.mut_type == 1.))
        nmito = int(torch.sum(self.mut_type == 2.))

        # create latent variables for parameter
        fpr,h_nucl,h_mito, u_nucl, lnu_nucl, u_mito,lnu_mito, u_cnv, lnu_cnv = torch.zeros(3), torch.zeros(1), torch.zeros(nmito), torch.zeros(nnucl), torch.zeros(nnucl), torch.zeros(1), torch.zeros(1), torch.tensor(1), torch.tensor(1)
        
        # hcnv and r can be a unique variable per chromosome or be celltype-specific
        if cnv_celltype:
            logith = torch.zeros(len(celltype_names)).repeat(ncnv,1)
            hcnv = torch.zeros(len(celltype_names)).repeat(ncnv,1)
            r = torch.zeros(len(celltype_names)).repeat(ncnv,1)
            cnv_ratio = torch.zeros(ncnv)
        else:
            hcnv = torch.zeros(nmut)
            r = torch.zeros(nmut)
        
        # VAF and CF 
        # allele frequencies can be present for different celltypes/timepoints and therefore are treated as different observations
        if class_af:
            
            cf_list = []
            
            # iterate through the different classes sample cf and compute VAF and fitted into a Binomial distribution
            for c in range(len(af_alpha)):
                
                # get class name
                nc = class_names[c]
                
                # sample clonal fraction
                cf_group = pyro.sample("cf_{}".format(nc), dist.Dirichlet(torch.ones(A.shape[0])))
                cf_list.append(cf_group)
                
                # compute theoretical allele frequencies
                af = get_af(cf_group, A)
                exome_counts = torch.round(af_alpha[c] + af_beta[c])
                
                for j in range(nmut):
                    
                    # get mutation name
                    m = names[j]
                    
                    # exome data is created by sampling with that allele frequency
                    pyro.sample("exome_{}_{}".format(nc, m), dist.Binomial(exome_counts[j], af[j]), obs = af_alpha[c][j])
                    
            # the final clonal fraction that is input it in the node assignmnet Categorical distribution is the weighted cf based on the abundances of 
            # different groups (e.g. 80% cancer and 20% T cells).
            ratio = torch.bincount(class_assign)/torch.tensor(class_assign.shape)
            cf = torch.sum(torch.stack(cf_list)*ratio.unsqueeze(1), dim = 0)
                
            
        # cases when there is only one sample and one measurement of bulk AFs
        else:
            # sample clonal fraction
            cf = pyro.sample("cf", dist.Dirichlet(torch.ones(A.shape[0])))

            # compute theoretical allele frequencies
            af = get_af(cf, A)
            exome_counts = torch.round(af_alpha + af_beta)
            
            for j in range(nmut):
                
                # get mutation name
                m = names[j]
                
                #exome data is created by sampling with that allele frequency
                pyro.sample("exome_{}".format(m), dist.Binomial(exome_counts[j], af[j]), obs = af_alpha[j])

        # FPR on for nuclear SNVs and one for mitochondrial SNVs
        for i in set(mut_type):

            i = int(i)

            n = types[i]
            
            # FPR parameter different for nuclear and mitochondrial mutations
            if i > 0:
                #fpr[i] = pyro.sample("fpr_{}".format(n), dist.Gamma(torch.tensor(40.), torch.tensor(5.)))
                #fpr[i] = pyro.sample("fpr_{}".format(n), dist.Gamma(torch.tensor(0.5), torch.tensor(5.)))
                fpr[i] = pyro.sample("fpr_{}".format(n), dist.Gamma(torch.tensor(7.), torch.tensor(200.)))
                
            # set unique heteroplasmy parameter for nuclear mutations
            if i == 1:
                
                h_nucl = pyro.sample("het_{}".format(n), dist.Beta(h_alpha[1], h_beta[1]))
                
            # set concetration parameter for mitochondrial mutations       
            if i == 2:

                # non-informative prior for concentration parameter
                u_mito = pyro.sample("u_{}".format(n), dist.Normal(0, 3))

                # exponentiate the overdispersion parameter 
                lnu_mito = torch.exp(u_mito)
                
            # set concentration parameter for CNVs
            if i == 0:
                # non-informative prior for concentration parameter
                u_cnv = pyro.sample("u_{}".format(n), dist.Normal(0, 3))

                # exponentiate the overdispersion parameter 
                lnu_cnv = torch.exp(u_cnv)
                

        # create mutation-specific latent variables (h, u and r)
        nucl_ind, mito_ind, cnv_ind = torch.tensor(0),torch.tensor(0), torch.tensor(0)
        for j in range(nmut):

            # get mutation name
            m = names[j]

            # add concentration parameters for nuclar variants
            if mut_type[j] == 1:
                               
                # concentration parameter for each nuclear variant
                u_nucl[nucl_ind] = pyro.sample("u_{}".format(m), dist.Normal(0, 3))
                
                lnu_nucl[nucl_ind] = torch.exp(u_nucl[nucl_ind])
                
                nucl_ind = torch.add(nucl_ind, torch.tensor(1))
                
            # add mutation-specific heteroplasmy for mito variants
            elif mut_type[j] == 2:
                
                h_mito[mito_ind] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha[2], h_beta[2]))
                
                mito_ind = torch.add(mito_ind, torch.tensor(1))
            
            # CNVs have an additional parameter r which represent the proportion of reads in the affected chr in cancer cells
            else:
                                            
                # celltype-specific CNV model
                if cnv_celltype:
                    for ind,ct in enumerate(celltype_names):
                                                
                        # latent variable logit h using celltype-specific priors inferred from patients w/o CNV
                        logith[cnv_ind][ind] = pyro.sample("hcnv_{}_{}".format(m, ct), dist.Normal(cnv_ct_mean[cnv_ind][ind], cnv_ct_sd[cnv_ind]))
                                          
                        # sigmoid transformation to compute hcnv (inverse of logit)
                        hcnv[cnv_ind][ind] = torch.sigmoid(logith[cnv_ind][ind])
                                                
                    # cnv_ratio in this case is the ratio between the fraction of reads in the affected chromosome between cancer and healthy
                    cnv_ratio[cnv_ind] = pyro.sample("cnv_ratio_{}".format(m), dist.Normal(r_cnv[j], 0.05))
               
                    # to compute the celltype-specific r hcvn is multiplied by r_ratio
                    #r[cnv_ind] = torch.mul(hcnv[cnv_ind], cnv_ratio)
                    
                    cnv_ind = torch.add(cnv_ind, torch.tensor(1))
                    
                # general model which assumes that the ratio of reads in the affected chromosome is the same in all cell populations        
                else:
                    # heteroplasmy as Beta distribution (high values e.g 1000 for nuclear SNVs, smaller for mitochondrial)
                    hcnv[j] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha[0], h_beta[0]))

                    # r in this case comes from a Beta distribution
                    #r[j] = pyro.sample("r_{}".format(m), dist.Beta(1, 1))
                    r[j] = pyro.sample("r_{}".format(m), dist.Beta(h_alpha[0], h_beta[0])) 

                    #the observation we have for this quantity is the ratio between r[j] (cancer) and hcnv[j] (healthy)
                    pyro.sample("cnv_state_{}".format(m), dist.Normal(torch.div(r[j],hcnv[j]).item(), 0.05), obs = r_cnv[j])
                    
        # to compute r we multiply the CNV-specific ratio and the healthy ratio hcnv
        if cnv_celltype:            
                    
            r = torch.mul(hcnv, cnv_ratio.unsqueeze(1))
                    
        # loop over all cells
        with pyro.plate("cell_loop", ncells):

            # marginalise node attachment using parallel enumeration
            # this means that for each cell it iterates over all attachments in parallel. 
            
            # whether to use bulk exome data or not
            if bulk_af:
                
                node = pyro.sample("assignment", dist.Categorical(probs = cf))
                
            else:
                
                node = pyro.sample("assignment", dist.Categorical(probs = torch.ones(A.shape[0])/A.shape[0]))

            # iterate through the mutations
            nucl_ind, mito_ind, cnv_ind = torch.tensor(0),torch.tensor(0), torch.tensor(0)
            for j in range(nmut):

                # get mutation name
                m = names[j]

                # transform data into tensor
                M_counts = M[:,j]
                N_counts = N[:,j]

                # make a tensor of booleans 
                k = torch.flatten(torch.gt(torch.where(A[node,j] == 1,1,0),0))
                
                # get mutation type (0: CNV, 1: nuclear, 2: mitochondrial)
                t = int(mut_type[j])

                # nuclear mutations
                if t == 1:             
                    # sample from the corresponding distribution 
                    pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                       dist.Poisson(fpr[t]), 
                                                                       dist.BetaBinomial(h_nucl*lnu_nucl[nucl_ind], (1-h_nucl)*lnu_nucl[nucl_ind], total_count = M_counts+N_counts)), 
                                                                        obs = M_counts)
                    
                    nucl_ind = torch.add(nucl_ind, torch.tensor(1))
                    
                # nuclear mutations
                elif t == 2:             
                    # sample from the corresponding distribution 
                    pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                       dist.Poisson(fpr[t]), 
                                                                       dist.BetaBinomial(h_mito[mito_ind]*lnu_mito, (1-h_mito[mito_ind])*lnu_mito, total_count = M_counts+N_counts)), 
                                                                        obs = M_counts)
                    
                    mito_ind = torch.add(mito_ind, torch.tensor(1))
                    
                # CNAs
                else:
                    if cnv_celltype:
                        
                        h_ct = torch.gather(hcnv[cnv_ind], 0, celltype)
                        r_ct = torch.clamp(torch.gather(r[cnv_ind], 0, celltype), 1e-10, 0.9999999999)
                                                              
                        # sample from beta binomial distribution
                        pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                           dist.BetaBinomial(h_ct*lnu_cnv, (1-h_ct)*lnu_cnv, total_count = N_counts),
                                                                           dist.BetaBinomial(r_ct*lnu_cnv, (1-r_ct)*lnu_cnv, total_count = N_counts)),
                                                                           obs = M_counts) 
                        cnv_ind = torch.add(cnv_ind, torch.tensor(1))
                    
                    else:
                        # sample from beta binomial distribution
                        pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                           dist.BetaBinomial(hcnv[j]*lnu_cnv, (1-hcnv[j])*lnu_cnv, total_count = N_counts),
                                                                           dist.BetaBinomial(r[j]*lnu_cnv, (1-r[j])*lnu_cnv, total_count = N_counts)),
                                                                           obs = M_counts) 

    # define guide
    def guide(self,M, N, A, mut_type, h_alpha, h_beta, af_alpha, af_beta, names, r_cnv, 
              bulk_af = True, class_af = False, class_names = [], class_assign = [],
              cnv_celltype = False, celltype = [], celltype_names = [], cnv_ct_mean = [], cnv_ct_sd = []):

        types = ["CNV", "nuclSNV", "mitoSNV"]

        nmut = M.shape[1]
        ncells = M.shape[0]
        nodes = A.shape[0]

        # get number of nuclear and mito mutations
        ncnv = int(torch.sum(self.mut_type == 0.))
        nnucl = int(torch.sum(self.mut_type == 1.))
        nmito = int(torch.sum(self.mut_type == 2.))
        
        # create latent variables for each mutation
        fpr,h_nucl, h_mito, u_nucl, u_mito, u_cnv = torch.zeros(3), torch.zeros(1), torch.zeros(nmito), torch.zeros(nnucl), torch.zeros(1), torch.zeros(1)
        
        # hcnv and r can be a unique variable per chromosome or be celltype-specific
        if cnv_celltype:
            logith = torch.zeros(len(celltype_names)).repeat(ncnv,1)
            hcnv = torch.zeros(len(celltype_names)).repeat(ncnv,1)
            cnv_ratio = torch.zeros(ncnv)
        else:
            h_cnv = torch.zeros(nmut)
            r = torch.zeros(nmut)

        # shape parameters of Dirichlet distribution for clonal fraction latent variable
        if class_af:
                        
            # iterate through the different classes and define parameters for cf Dirichlet distribution
            for c in range(len(af_alpha)):
                
                # get class name
                nc = class_names[c]
                
                # parameters of Dirichlet distribution for particular class
                cf_param = pyro.param("cf_conc_{}".format(nc), torch.ones(A.shape[0]), constraint = constraints.positive)
                
                # cf latent variable for specific class
                cf = pyro.sample("cf_{}".format(nc), dist.Dirichlet(cf_param))

        # cases when there is only one sample and one measurement of bulk AFs
        else:
            # define parameters of Dirichlet distribution
            cf_param = pyro.param("cf_conc", torch.ones(A.shape[0]), constraint = constraints.positive)
            cf = pyro.sample("cf", dist.Dirichlet(cf_param))

        for i in set(mut_type):

            i = int(i)

            n = types[i]
            
            # set hyperparameters of nuclear heteroplasmy parameter
            if i == 1:
                
                # mean and concentration of the beta distribution as parameters
                h_mean_param = pyro.param("h_mean_{}".format(n), torch.tensor(0.5))
                h_conc_param = pyro.param("h_conc_{}".format(n), torch.tensor(3.0))

                # transform to shape parameters (helps to avoid non-identifiabilities)
                h_alpha = torch.exp(h_conc_param) * torch.sigmoid(h_mean_param)
                h_beta = torch.exp(h_conc_param) * (1-torch.sigmoid(h_mean_param))

                # beta distribution to approximate posterior of heteroplasmy latent variable
                h_nucl = pyro.sample("het_{}".format(n), dist.Beta(h_alpha, h_beta))
                
            # set hyperparameters of mito concentration parameters
            if i == 2:
                
                # mean and sd of Normal distribution for concentration parameters
                u_mean_param = pyro.param("u_mean_{}".format(n), torch.tensor(0.))
                u_sd_param = pyro.param("u_sd_{}".format(n), torch.tensor(3.0), constraint = constraints.positive)                      

                # normal distribution to approximate concentration parameters
                u_mito = pyro.sample("u_{}".format(n), dist.Normal(u_mean_param, u_sd_param))
                
            # set hyperparameters of CNV concentration parameter
            if i == 0:
                
                # mean and sd of Normal distribution for concentration parameters
                u_mean_param = pyro.param("u_mean_{}".format(n), torch.tensor(0.))
                u_sd_param = pyro.param("u_sd_{}".format(n), torch.tensor(3.0), constraint = constraints.positive)                      

                # normal distribution to approximate concentration parameters
                u_cnv = pyro.sample("u_{}".format(n), dist.Normal(u_mean_param, u_sd_param))

            # set hyperparameters of mito and nuclear FPR
            if i > 0:
                fpr_shape = pyro.param("fpr_shape_{}".format(n), torch.log(torch.tensor(5)), constraint = constraints.positive)
                fpr_scale_param = pyro.param("fpr_scale_{}".format(n), torch.tensor(0.2), constraint = constraints.positive)
                fpr_scale = torch.div(1, fpr_scale_param)

                # gamma dist. to approximate posterior of fpr latent variable
                fpr = pyro.sample("fpr_{}".format(n), dist.Gamma(torch.exp(fpr_shape), fpr_scale))

        # loop over the different mutation types
        for j in range(nmut):

            # get mutation name
            m = names[j]
            
            nucl_ind, mito_ind, cnv_ind = torch.tensor(0),torch.tensor(0),torch.tensor(0)
            
            # set hyperparameters of nuclear mutation concentration parameter
            if mut_type[j] == 1:
                
                # mean and sd of Normal distribution for concentration parameters
                u_mean_param = pyro.param("u_mean_{}".format(m), torch.tensor(0.))
                u_sd_param = pyro.param("u_sd_{}".format(m), torch.tensor(3.0), constraint = constraints.positive)                      

                # normal distribution to approximate concentration parameters
                u_nucl[nucl_ind] = pyro.sample("u_{}".format(m), dist.Normal(u_mean_param, u_sd_param))
                
                nucl_ind = torch.add(nucl_ind, torch.tensor(1))
                
            # set hyperparameter of mito mutations heteroplasmy
            if mut_type[j] == 2:
                
                # mean and concentration of the beta distribution as parameters
                h_mean_param = pyro.param("h_mean_{}".format(m), torch.tensor(0.5))
                h_conc_param = pyro.param("h_conc_{}".format(m), torch.tensor(3.0))

                # transform to shape parameters (helps to avoid non-identifiabilities)
                h_alpha = torch.exp(h_conc_param) * torch.sigmoid(h_mean_param)
                h_beta = torch.exp(h_conc_param) * (1-torch.sigmoid(h_mean_param))

                # beta distribution to approximate posterior of heteroplasmy latent variable
                h_mito[mito_ind] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha, h_beta))
                
                mito_ind = torch.add(mito_ind, torch.tensor(1))
                
                   
            # for CNVs parameters for logitr and logithcnv are created
            if mut_type[j] == 0:
                if cnv_celltype:
                    for ind, ct in enumerate(celltype_names):
                        
                        # mean and sd for celltype specific hcnv parameters
                        hcnv_mean_param = pyro.param("hcnv_mean_{}.{}".format(m, ct), cnv_ct_mean[cnv_ind][ind])
                        hcnv_sd_param = pyro.param("hcnv_sd_{}.{}".format(m, ct), cnv_ct_sd[cnv_ind], constraint = constraints.positive)
                        
                        # latent variable hcnv
                        logith[cnv_ind][ind] = pyro.sample("hcnv_{}_{}".format(m, ct), dist.Normal(hcnv_mean_param, hcnv_sd_param))
                        
                    # mean and sd of Normal distribution of cnv_ratio latent variable
                    cnv_ratio_mean = pyro.param("CNVratio_mean_{}".format(m), r_cnv[j])
                    cnv_ratio_sd = pyro.param("CNVratio_sd_{}".format(m), torch.tensor(0.05), constraint = constraints.positive)

                    # latent variable cnv_ratio
                    cnv_ratio[cnv_ind] = pyro.sample("cnv_ratio_{}".format(m), dist.Normal(cnv_ratio_mean, cnv_ratio_sd))
                    
                    cnv_ind = torch.add(cnv_ind, torch.tensor(1))
                        
                                            
                else:
                    
                    # mean and concentration of the beta distribution as parameters
                    h_mean_param = pyro.param("h_mean_{}".format(m), torch.tensor(0.05), constraint = constraints.unit_interval)
                    h_conc_param = pyro.param("h_conc_{}".format(m), torch.tensor(2.), constraint = constraints.positive)

#                     # transform to shape parameters (helps to avoid non-identifiabilities)
#                     h_alpha = torch.exp(h_conc_param) * torch.sigmoid(h_mean_param)
#                     h_beta = torch.exp(h_conc_param) * (1-torch.sigmoid(h_mean_param))
                    # transform parameters 
                    h_alpha = torch.mul(h_mean_param,h_conc_param)
                    h_beta = torch.mul((1-h_mean_param),h_conc_param)

                    # beta distribution to approximate posterior of heteroplasmy latent variable
                    h_cnv[j] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha, h_beta))
                    
                    # mean and concentration of the beta distribution as parameters
                    r_mean_param = pyro.param("r_mean_{}".format(m), torch.tensor(0.05), constraint = constraints.unit_interval)
                    r_conc_param = pyro.param("r_conc_{}".format(m), torch.tensor(2.), constraint = constraints.positive)

                    # beta distribution to approximate posterior of heteroplasmy latent variable
#                     r_alpha = torch.exp(r_conc_param) * torch.sigmoid(r_mean_param)
#                     r_beta = torch.exp(r_conc_param) * (1-torch.sigmoid(r_mean_param))

                    r_alpha = torch.mul(r_mean_param,r_conc_param)
                    r_beta = torch.mul((1-r_mean_param),r_conc_param)

                    # latent variable r
                    r[j] = pyro.sample("r_{}".format(m), dist.Beta(r_alpha, r_beta))  
                
    # define guide for clonal assignment
    @config_enumerate
    def guide_clones(self, M, N, A, mut_type, h_alpha, h_beta, af_alpha, af_beta, names, r_cnv, 
                     bulk_af = True, class_af = False, class_names = [], class_assign = [],
                     cnv_celltype = False, celltype = [], celltype_names = [], cnv_ct_mean = [], cnv_ct_sd = []):
        
        # get params to block 
        params = [i for i in get_params(1)['svi_params'].keys()]

        # get number of cells and nodes
        ncells = M.shape[0]
        nodes = A.shape[0]

        # use the inferred parameters for the tree guide
        with poutine.block(hide = params):
            self.guide(M = M,
                     N = N,
                     mut_type = mut_type,
                     A = A,
                     h_alpha = h_alpha,
                     h_beta = h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = r_cnv,
                     names = names,
                     class_af = class_af,
                     bulk_af = bulk_af,
                     class_names = class_names,
                     class_assign = class_assign, 
                     cnv_celltype = cnv_celltype, 
                     celltype = celltype,
                     celltype_names = celltype_names, 
                     cnv_ct_mean = cnv_ct_mean, 
                     cnv_ct_sd = cnv_ct_sd)

        # parameter for clonal assignment
        with pyro.plate("cell_loop", ncells):

            # parameter for clonal assignment
            clone_probs = pyro.param("clonal_probs", torch.ones(ncells, nodes)/nodes, constraint=constraints.simplex)
            
            # latent variable for assigned node
            pyro.sample("assignment", dist.Categorical(clone_probs))  
            
                 
    # function to select best trees after adding a new mutation
    def select_tree(self, num_iter, init, num_particles = 5, vect_particles=False):
        
        print("Choosing best tree among {} potential candidates".format(len(self.potential_trees)))
        
        start = timer()
        
        self.num_iter = num_iter

        optimiser = pyro.optim.AdagradRMSProp({})
        # loss function
        loss_func = pyro.infer.TraceEnum_ELBO(max_plate_nesting = 1, num_particles = num_particles, vectorize_particles = vect_particles)

        # set up inference algorithm
        svi = SVI(self.model, self.guide, optimiser, loss=loss_func)
        
        # create dictionaries to store values
        elbo = {}
        step_params = {}
        best_iter = {}
        
        # dictionary to count how many times each tree is selected    
        tree_count = torch.zeros(len(self.potential_trees))
        
        # subsetting allele frequencies is different if there are several populations or not
        if(self.class_af):
            
            af_alpha = self.af_alpha[:,self.muts]
            af_beta = self.af_beta[:,self.muts]
            
        else:
            af_alpha = self.af_alpha[self.muts]
            af_beta = self.af_beta[self.muts]
        
        # run the model for all potential trees
        for tree in range(len(self.potential_trees)):
            
            tree_start = timer()

            # clear previous inferred parameters
            pyro.clear_param_store()
            elbo[tree] = torch.zeros(num_iter)
            step_params[str(tree)] = {}
            
            for i in range(num_iter):
               
                elbo[tree][i] = svi.step(M = self.M[:, self.muts],
                                         N = self.N[:, self.muts],
                                         mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                                         A = self.potential_trees[tree],
                                         h_alpha = self.h_alpha,
                                         h_beta = self.h_beta,
                                         af_alpha = af_alpha,
                                         af_beta = af_beta,
                                         r_cnv = self.r_cnv[self.muts],
                                         names = [self.names[i] for i in self.muts],
                                         bulk_af = self.bulk_af,
                                         class_af = self.class_af,
                                         class_names = self.class_names, 
                                         class_assign = self.class_assign,
                                         cnv_celltype = self.cnv_celltype,
                                         celltype = self.celltype,
                                         celltype_names = self.celltype_names,
                                         cnv_ct_mean = self.cnv_ct_mean,
                                         cnv_ct_sd = self.cnv_ct_sd)                
                # store parameter values
                params = get_params(1)['svi_params']
                for key, value in params.items():
                    if i == 0:
                        step_params[str(tree)][str(key)] = []
                    step_params[str(tree)][str(key)].append(value)
                    
                # remove previously computed posterior predictive distributions
                if tree == 0:
                    self.post_predictive = {}              
                    
            # store iteration with lowest ELBO
            best_iter[tree] = torch.argmin(elbo[tree])
            
            # print run time
            tree_end = timer()
            tree_time = timedelta(seconds=tree_end-tree_start)
            m, s = divmod(tree_time.seconds, 60)
            print("Model fit to tree {} in {}m {}s".format(tree, m, s))

        # determine trees with lowest ELBO by taking any tree which has the lowest ELBO in at least one of the last 50 iterations
        elbos = torch.zeros(len(self.potential_trees))
        best_trees = []
        for it in range(50):
            elbos = torch.zeros(len(self.potential_trees))
            for t in range(len(self.potential_trees)):
                elbos[t] = elbo[t][num_iter-it-1]
            best_trees.append(torch.argmin(elbos).item())
                
                          
        # add best tree(s) as attribute for further building
        self.tree = [self.potential_trees[i] for i in list(set(best_trees))]
        self.tree_indices = list(set(best_trees))
        self.best_iter = best_iter
        
        # print output
        if len(self.tree) == 1:
            print("Tree {} selected for mutations {}:".format(self.tree_indices, [self.names[i] for i in self.muts]))
            print(self.tree)
        else:
            print("Trees {} selected for mutations {}".format(self.tree_indices, [self.names[i] for i in self.muts]))
            for tree in self.tree:
                print(tree)
        
        # print run time
        end = timer()
        tim = timedelta(seconds=end-start)
        h, r = divmod(tim.seconds, 3600)
        m, s = divmod(tim.seconds, 60)
        print("Total run time {}h {}m {}s".format(h, m, s))
        
        # add elbo and parameters as attribute
        self.elbo = elbo
        self.parameters = step_params
        
        # print line plot with ELBOs
        self.print_elbo(num_iter = num_iter, init = init)
        
    # function to run the model on a single tree
    def single_tree(self, tree, mut_names, num_iter = 0):
        
        # run svi for the tree of interest        
        optimiser = pyro.optim.AdagradRMSProp({})
                  
        # compute clonal assignments
        loss_func = pyro.infer.TraceEnum_ELBO(max_plate_nesting = 1, num_particles = 5)

        # Infer parameters for the tree of interest
        svi = SVI(self.model, self.guide, optimiser, loss=loss_func)
        
        # clear previous parameters
        pyro.clear_param_store()
        
        # Generate children and parents from the current tree
        for i in tree.shape[1]-1:
            
            if i == 0:
                break
            else:
                self.add_mutation()
                
        # get index of the selected tree
        tree_ind = self.potential_trees.index(tree)
        
        # save children and parent options
        self.children = self.children[tree_ind]
        self.parents = self.parents[tree_ind]
        
        # add mutation names
        self.muts = mut_names
        
        # add tree as selected one
        self.tree, self.potential_trees = tree
        self.tree_indices = 0
        
        if num_iter == 0:
            num_iter = self.num_iter

        # iterate until lowest ELBO was obtained
        init = self.best_iter[tree]
        
        # subsetting allele frequencies is different if there are several populations or not
        if(self.class_af):
            
            af_alpha = self.af_alpha[:,self.muts]
            af_beta = self.af_beta[:,self.muts]
            
        else:
            af_alpha = self.af_alpha[self.muts]
            af_beta = self.af_beta[self.muts]
            
        # run model with clonal assignments marginalised out    
        for i in range(init):
            
            svi.step(M = self.M[:, self.muts],
                     N = self.N[:, self.muts],
                     mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                     A = tree,
                     h_alpha = self.h_alpha,
                     h_beta = self.h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
                     bulk_af = self.bulk_af,
                     class_af = self.class_af,
                     class_names = self.class_names, 
                     class_assign = self.class_assign,
                     cnv_celltype = self.cnv_celltype,
                     celltype = self.celltype,
                     celltype_names = self.celltype_names,
                     cnv_ct_mean = self.cnv_ct_mean,
                     cnv_ct_sd = self.cnv_ct_sd)
        
        print("Computing clonal assignment probabilities for tree {}".format(tree))

        # set up inference algorithm with guide to compute clonal probabilities
        svi = SVI(self.model, self.guide_clones, optimiser, loss=loss_func)
        
        start = timer()
        
        # run model to compute clonal probabilities (fixing the parameters learnt in the former step)
        for i in range(num_iter):
            
            svi.step(M = self.M[:, self.muts],
                     N = self.N[:, self.muts],
                     mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                     A = self.tree[tree],
                     h_alpha = self.h_alpha,
                     h_beta = self.h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
                     bulk_af = self.bulk_af,
                     class_af = self.class_af,
                     class_names = self.class_names, 
                     class_assign = self.class_assign,
                     cnv_celltype = self.cnv_celltype,
                     celltype = self.celltype,
                     celltype_names = self.celltype_names,
                     cnv_ct_mean = self.cnv_ct_mean,
                     cnv_ct_sd = self.cnv_ct_sd)
    
        # print run time
        tree_end = timer()
        tree_time = timedelta(seconds=tree_end-start)
        m, s = divmod(tree_time.seconds, 60)
        
        print("Clonal assignment probabilities computed in {}m and {}s".format(m,s))
        
        # save clonal probabilites as attribute
        if len(self.clone_probs) == 0:
            
            self.clone_probs = {}
            
        self.clone_probs[1] = pyro.param("clonal_probs")
        
                
    # function to add mutation in all possible sites of the tree
    def add_mutation(self):

        new_trees = []
        children_list = []
        parents_list = []

        for i in range(len(self.potential_trees)):
            
            if i not in self.tree_indices:
                continue

            tree = self.potential_trees[i]

            nmuts = tree.shape[1]
            nodes = tree.shape[0]

            # iterate through the nodes
            for node in range(nodes):

                # add mutation column to the current tree
                new_tree = torch.cat((tree, torch.zeros(tree.shape[0],1)), dim =1)

                # if node is root create tree with mutation indpendent from all others or upstream of all others
                if node == 0:

                    # Add mutation as independent node hanging from the root -------------------------------------
                    # create new node with only the mutation of interest
                    new_node = torch.zeros(nmuts+1)
                    new_node[nmuts] = 1

                    new_trees.append(torch.cat((new_tree, new_node.unsqueeze(dim=0)), dim = 0))

                    # create entry for children nodes for new tree
                    child = copy.deepcopy(self.children[i])
                    child[0].append(nodes)
                    child.append([])
                    children_list.append(child)

                    # create entry for parent nodes in the new tree
                    parent = copy.deepcopy(self.parents[i])
                    parent.append([0])
                    parents_list.append(parent)

                    # add mutation upstream of the rest -----------------------------------------------------
                    # add mutation column to the current tree
                    new_tree = torch.cat((tree, torch.zeros(tree.shape[0],1)), dim =1)

                    # create new node
                    new_node = torch.zeros(nmuts+1)
                    new_node[nmuts] = 1

                    # add mutation to all nodes except the root
                    new_tree[1:,nmuts] = 1

                    # add new tree to list
                    new_trees.append(torch.cat((new_tree, new_node.unsqueeze(dim=0)), dim = 0))

                    # create entry in parent and children nodes for new tree
                    aff_nodes = [i for i in range(nodes)[1:]]
                    child = copy.deepcopy(self.children[i])
                    pt = copy.deepcopy(self.parents[i])
                    child[0].append(nodes)

                    for p in aff_nodes:
                        pt[p].append(nodes)

                    child.append(aff_nodes)
                    pt.append([0])
                    children_list.append(child)
                    parents_list.append(pt)


                # for all other nodes I create a tree with the mutation upstream, merged with the current node and downstream.
                else:
                    
                    # add mutation upstream of current node ----------------------------------------------
                    # add mutation column to the current tree
                    new_tree = torch.cat((tree, torch.zeros(tree.shape[0],1)), dim =1)

                    # get parent nodes 
                    pt = self.parents[i][node] 

                    # get immediate parent node (node with most mutations among the parents)
                    pt_im = tree[pt,:][torch.argsort(torch.sum(tree[pt,:], dim = 1)), :][-1,:] 

                    # create new node adding new mutation to immediate parent node
                    new_node = torch.cat((pt_im, torch.tensor(1).unsqueeze(dim = 0)), dim = 0)

                    # attach new node to tree
                    new_tree = torch.cat((new_tree, new_node.unsqueeze(dim = 0)), dim = 0)

                    # get children nodes 
                    ch = copy.deepcopy(self.children[i][node])
                    ch.append(node)

                    # add mutation to children nodes in the tree
                    new_tree[ch, nmuts] = 1

                    # check whether tree exists 
                    for tr in new_trees:
                        if torch.equal(new_tree,tr):
                            exists = True
                        else:
                            exists = False

                    if not exists:

                        # add new tree to list
                        new_trees.append(new_tree)

                        # create children list for tree
                        child = copy.deepcopy(self.children[i])
                        # add current node as children of parent nodes
                        for p in pt:

                            child[p].append(nodes)

                        # add entry of children nodes
                        child.append(ch)

                        children_list.append(child)

                        # create parent list for tree
                        parent = copy.deepcopy(self.parents[i])
                        for p in ch:

                            parent[p].append(nodes)

                        # add entry for parents of newly added node
                        parent.append(pt)

                        # add parents list
                        parents_list.append(parent)

                    # merge mutation to the current node ----------------------------------------------------------
                    # add mutation column to the current tree
                    new_tree = torch.cat((tree, torch.zeros(tree.shape[0],1)), dim =1)

                    # get children nodes 
                    ch = copy.deepcopy(self.children[i][node])

                    # add mutation to children nodes in the tree
                    new_tree[ch, nmuts] = 1

                    # add mutation merged to current node
                    new_tree[node, nmuts] = 1

                    # add new tree to list
                    new_trees.append(new_tree)

                    # create children list for tree. As we don't add any new node it stays the same        
                    children_list.append(copy.deepcopy(self.children[i]))

                    # create parents list for new tree. It also stays the same
                    parents_list.append(copy.deepcopy(self.parents[i]))

                    # add mutation as terminal node -------------------------------------------------------------------
                    # add mutation column to the current tree
                    new_tree = torch.cat((tree, torch.zeros(tree.shape[0],1)), dim =1)

                    # get current node
                    curr_node = copy.deepcopy(tree[node,:])

                    # add mutation downstream of current node
                    new_node = torch.cat((curr_node, torch.tensor(1).unsqueeze(dim = 0)), dim = 0)

                    # attach new node to tree
                    new_tree = torch.cat((new_tree, new_node.unsqueeze(dim = 0)), dim = 0)  

                    # add new tree to list
                    new_trees.append(new_tree)

                    # get parents of current node
                    pt = copy.deepcopy(self.parents[i][node])
                    pt.append(node)

                    # create children list for new tree
                    child = copy.deepcopy(self.children[i])
                    for p in pt:

                        child[p].append(nodes)

                    # add entry for new terminal node with children
                    child.append([])

                    children_list.append(child)

                    # create parent list for new tree
                    parent = copy.deepcopy(self.parents[i])
                    parent.append(pt)

                    parents_list.append(parent)    

                    # For branching nodes an additional tree should be created with the new mutation downstream of the node
                    # but upstream of the mutually exclusive nodes.
                    if node < nodes-1:
                        
                        # check whether node is branching point by comparing number of children nodes with immediate downstream node
                        curre_ch = len(self.children[i][node])
                        dw_ch_list = copy.deepcopy(self.children[i][node+1:])
                        dw_ch_list.sort(key = len)
                        dw_ch = len(dw_ch_list[-1])
                        
                        
                        if curre_ch-dw_ch > 1:

                            # add mutation column to the current tree
                            new_tree = torch.cat((tree, torch.zeros(tree.shape[0],1)), dim =1)

                            # get current node
                            curr_node = copy.deepcopy(tree[node,:])

                            # add mutation downstream of current node
                            new_node = torch.cat((curr_node, torch.tensor(1).unsqueeze(dim = 0)), dim = 0)

                            # attach new node to tree
                            new_tree = torch.cat((new_tree, new_node.unsqueeze(dim = 0)), dim = 0)  

                            # get children nodes
                            ch = copy.deepcopy(self.children[i][node])

                            # add mutation to children nodes
                            new_tree[ch, nmuts] = 1

                            # add new tree to list
                            new_trees.append(new_tree)

                            # create children list for tree
                            child = copy.deepcopy(self.children[i])
                            pt = copy.deepcopy(self.parents[i][node])
                            # add current node as children of parent nodes
                            for p in pt:

                                child[p].append(nodes)

                            # add entry of children nodes
                            child.append(ch)

                            children_list.append(child)

                            # create parent list for tree
                            parent = copy.deepcopy(self.parents[i])
                            for p in ch:

                                parent[p].append(nodes)

                            pt.append(node)
                            parent.append(pt)

                            # add parents list
                            parents_list.append(parent)


        # make sure no duplicated trees are present
        incl_ind = []
        tree_list = []
        for ind, i in enumerate(new_trees):
            
            if any((i== t_).all() for t_ in tree_list if i.shape == t_.shape):
                continue
            else:    
                tree_list.append(i)
                incl_ind.append(ind)
            
        self.potential_trees = [new_trees[i] for i in incl_ind]
        self.children = [children_list[i] for i in incl_ind]
        self.parents = [parents_list[i] for i in incl_ind]
                         
        # choose mutation to add based on coverage
        mut = self.sel_mutation()
        
        print("{} added to the tree".format(self.names[mut]))

    
    # function to plot the ELBOs of trees during SVI
    def print_elbo(self, num_iter, init, include = []):

        data = pd.DataFrame(data = {"ELBO": [i.item() for key in self.elbo for i in self.elbo[key][init:num_iter]],
                                    "tree": itertools.chain.from_iterable(zip(*itertools.repeat(range(len(self.potential_trees)), num_iter-init))),
                                    "num_iter": [*range(num_iter-init)]*len(self.potential_trees)})
        # reduce high values
        data.assign(ELBO = [1.5e5 if ELBO > 1.5e5 else ELBO for ELBO in data['ELBO']])
        
        if len(include) > 0:
            sns.lineplot(data = data[data["tree"].isin(include)], x = "num_iter", y = "ELBO",hue="tree", palette = sns.color_palette("tab10", len(include)))
            
        else:
            sns.lineplot(data = data, x = "num_iter", y = "ELBO", hue="tree", palette = sns.color_palette("tab10", len(self.potential_trees)))
    
    # function to print the parameter values over the SVI iterations 
    # it indicates the number of initial iterations to exclude
    def print_params(self, include = [], it = 0):
               
        params = self.parameters["0"].keys()

        sel_trees = include

        data_params = pd.DataFrame({"tree": itertools.chain.from_iterable(zip(*itertools.repeat(sel_trees, self.num_iter))),
                                    "num_iter": [*range(self.num_iter)]*len(sel_trees)})
        for p in params:

            if p.startswith("cf_conc"):

                data_params[p+"_0"] = [i[0] for key in sel_trees for i in self.parameters[str(key)][p]]

            else:
                data_params[p] = [i.item() for key in sel_trees for i in self.parameters[str(key)][p]]

        rows = math.ceil(len(params)/4)+1

        fig, ax = plt.subplots(nrows = rows, ncols = 4, figsize = (35, math.ceil(rows)*6))
        count = 0
        
        if it == 0:
            for index, p in enumerate(params):

                count += 1
                row = divmod(count, 4)[0]
                col = divmod(count, 4)[1]

                if p.startswith("cf_conc"):
                    
                    name = p+ "_0"
                    
                    sns.lineplot(data = data_params, x = "num_iter", y = name, hue="tree", palette = sns.color_palette("tab10", len(sel_trees)), ax=ax[row, col])
                    ax[row, col].set(title = name, xlabel = "iterations")

                else:
                    sns.lineplot(data = data_params, x = "num_iter", y = p, hue="tree", palette = sns.color_palette("tab10", len(sel_trees)), ax=ax[row, col])
                    ax[row, col].set(title = p, xlabel = "iterations")
                    
        else:
            for index, p in enumerate(params):

                count += 1
                row = divmod(count, 4)[0]
                col = divmod(count, 4)[1]
                
                if p.startswith("cf_conc"):
                    
                    name = p + "_0"
                    
                    sns.lineplot(data = data_params[data_params["num_iter"] > it], x = "num_iter", y = name, hue="tree", 
                                 palette = sns.color_palette("tab10", len(sel_trees)), ax=ax[row, col])
                    ax[row, col].set(title = name, xlabel = "iterations")

                else:
                    sns.lineplot(data = data_params[data_params["num_iter"] > it], x = "num_iter", y = p, hue="tree", 
                                 palette = sns.color_palette("tab10", len(sel_trees)), ax=ax[row, col])
                    ax[row, col].set(title = p, xlabel = "iterations")            

    # function to retrieve the parameters of the posterior distributions
    def get_post_params(self):
        
        params = {}
        
        for ind,tr in enumerate(self.tree_indices):
            
            params[tr] = {}
            
            best_iter = self.best_iter[tr]
            
            for p in self.parameters['0'].keys():
                               
                params[tr][p] = self.parameters[str(tr)][p][best_iter]
                
                # compute clonal fractions and allele frequencies
                if p.startswith("cf_conc_"):
                    
                    cf_name = p.split("_")[0] + "_" + p.split("_")[2]
                    af_name = "af_" + p.split("_")[2]
                    
                    cf = torch.Tensor(params[tr][p])/torch.sum(torch.Tensor(params[tr][p]))
                    params[tr][cf_name] = cf.detach().cpu().numpy()
                    params[tr][af_name] = get_af(cf, self.tree[ind]).detach().cpu().numpy()
                
                # allele frequencies for sample with only one set of AF priors(no different populations/timepoints)
                if p.startswith("cf_conc_"):
                    
                    cf_name = p.split("_")[0] 
                    af_name = "af"
                    
                    cf = torch.Tensor(params[tr][p])/torch.sum(torch.Tensor(params[tr][p]))
                    params[tr][cf_name] = cf.detach().cpu().numpy()
                    params[tr][af_name] = get_af(cf, self.tree[ind]).detach().cpu().numpy()
            
            
        self.params_post = params
        return(params)
    
    # function to generate and approximate posterior distribution given learnt parameters
    def get_posterior(self, variable = [], tree = [], data_points = 1000):
                
        data = np.zeros(data_points)
        
        # simulate posterior from a Gamma a distribution
        if variable.startswith("fpr"):
        
            ty = variable.split("_")[1]
            
            for i in range(data_points):
                data[i] = dist.Gamma(np.exp(self.params_post[tree]["fpr_shape_" + ty].item()),
                                     1/(self.params_post[tree]["fpr_scale_" + ty].item())).sample()
        
        # normal distribution for the concentration parameter
        elif variable.startswith("u") or variable.startswith("CNVratio"):
            
            vname = variable.split("_")[0]
            ty = variable.split("_")[1]
            
            for i in range(data_points):
                data[i] = dist.Normal(self.params_post[tree][vname+"_mean_" + ty].item(),
                                      self.params_post[tree][vname+"_sd_" + ty].item()).sample()
                
        elif (variable.startswith("hcnv") or variable.startswith("r")) and self.cnv_celltype:
            
            vname = variable.split("_")[0]
            ty = variable.split("_")[1]
            
            for i in range(data_points):
                data[i] = dist.Normal(self.params_post[tree][vname + "_mean_" + ty].item(),
                                      self.params_post[tree][vname + "_sd_" + ty].item()).sample()
                
            data = torch.sigmoid(torch.Tensor(data)).detach().cpu().numpy()            
                
        # Beta distribution for h and r latent variables
        elif variable.startswith("h_chr") or variable.startswith("r"):
            
            v, mut = variable.split("_")
            
            
            mean = torch.tensor(self.params_post[tree][v + "_mean_" + mut].item())
            conc = torch.tensor(self.params_post[tree][v + "_conc_" + mut].item())
            
            h_alpha = torch.mul(mean,conc)
            h_beta = torch.mul((1-mean),conc)

            
            for i in range(data_points):
                data[i] = dist.Beta(h_alpha, h_beta).sample()
                
        # beta distribution for heteroplasmy latent variables but with different re-parameterization
        else:

            v, mut = variable.split("_")

            mean = torch.tensor(self.params_post[tree][v + "_mean_" + mut].item())
            conc = torch.tensor(self.params_post[tree][v + "_conc_" + mut].item())

            h_alpha = torch.exp(conc) * torch.sigmoid(mean)
            h_beta = torch.exp(conc) * (1-torch.sigmoid(mean))


            for i in range(data_points):
                data[i] = dist.Beta(h_alpha, h_beta).sample()
                
        return(data)
                    
    # function to plot the approximated posterior distribution of the latent variables for one tree
    # options are: fpr,h,r,u, or all
    def plot_posterior(self, variables = ["fpr_"], tree = [], data_points = 1000):
        
        # retrieve parameters of approximated posterior distributions
        self.get_post_params()
               
        # if all variables are selected then plot all variables
        if variables == "all":
            variables = ["fpr_", "h_", "u_", "r_", "hcnv_", "CNVratio_"]
        
        # iterate through selected variables
        for v in variables:
                      
            # get parameter names
            param_names = [i for i in self.parameters['0'].keys() if i.startswith(v)]
            
            # check if params found (if all are selected, some tree will not have r)
            if len(param_names) == 0:
                continue
            
            # get number of plots to show
            num_plots = int(len(param_names)/2)
            
            if(v == "hcnv" and self.cnv_celltype):
                
                row = int(torch.sum(self.mut_type == 0.))*2
                
            elif v == "hcnv":
                
                row = int(torch.sum(self.mut_type == 0.))
                
            else:
                
                row = 1
            
            fig, ax = plt.subplots(nrows = row, ncols = num_plots, figsize = 
                                   (int(num_plots*5), 5))
                
            for i in range(num_plots):

                name = v + param_names[i*2].split("_")[2]
                                    
                data = self.get_posterior(name, tree, data_points)
                
                if num_plots > 1:                
                    sns.histplot(data, ax = ax[i])
                    ax[i].set(title = name)
                else:
                    sns.histplot(data).set_title(name)
                                                    
    # function to compute posterior distribution of clonal assignment
    def clonal_assignment(self, tree, num_iter = 0):
        
        # run svi for the tree of interest        
        optimiser = pyro.optim.AdagradRMSProp({})
                  
        # compute clonal assignments
        loss_func = pyro.infer.TraceEnum_ELBO(max_plate_nesting = 1, num_particles = 5)

        # Infer parameters for the tree of interest
        svi = SVI(self.model, self.guide, optimiser, loss=loss_func)
        
        # clear previous parameters
        pyro.clear_param_store()
        
        if num_iter == 0:
            num_iter = self.num_iter

        # iterate until lowest ELBO was obtained
        nit = self.best_iter[tree]
        
        # subsetting allele frequencies is different if there are several populations or not
        if(self.class_af):
            
            af_alpha = self.af_alpha[:,self.muts]
            af_beta = self.af_beta[:,self.muts]
            
        else:
            af_alpha = self.af_alpha[self.muts]
            af_beta = self.af_beta[self.muts]
            
        for i in range(nit):
            
            svi.step(M = self.M[:, self.muts],
                     N = self.N[:, self.muts],
                     mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                     A = self.tree[tree],
                     h_alpha = self.h_alpha,
                     h_beta = self.h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
                     bulk_af = self.bulk_af,
                     class_af = self.class_af,
                     class_names = self.class_names, 
                     class_assign = self.class_assign,
                     cnv_celltype = self.cnv_celltype,
                     celltype = self.celltype,
                     celltype_names = self.celltype_names,
                     cnv_ct_mean = self.cnv_ct_mean,
                     cnv_ct_sd = self.cnv_ct_sd)
        
        print("Computing clonal assignment probabilities for tree {}".format(tree))

        # set up inference algorithm
        svi = SVI(self.model, self.guide_clones, optimiser, loss=loss_func)
        
        start = timer()
        
        for i in range(num_iter):
            
            svi.step(M = self.M[:, self.muts],
                     N = self.N[:, self.muts],
                     mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                     A = self.tree[tree],
                     h_alpha = self.h_alpha,
                     h_beta = self.h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
                     bulk_af = self.bulk_af,
                     class_af = self.class_af,
                     class_names = self.class_names, 
                     class_assign = self.class_assign,
                     cnv_celltype = self.cnv_celltype,
                     celltype = self.celltype,
                     celltype_names = self.celltype_names,
                     cnv_ct_mean = self.cnv_ct_mean,
                     cnv_ct_sd = self.cnv_ct_sd)
    
        # print run time
        tree_end = timer()
        tree_time = timedelta(seconds=tree_end-start)
        m, s = divmod(tree_time.seconds, 60)
        
        print("Clonal assignment probabilities computed in {}m and {}s".format(m,s))
        
        # save clonal probabilites as attribute
        if len(self.clone_probs) == 0:
            
            self.clone_probs = {}
            
        self.clone_probs[tree+1] = pyro.param("clonal_probs")
        
        
    # function to export potential_trees, parameters of posterior distribution and clonal assignment probabilities
    def export_pickle(self, file):
        
        obj = {}
        
        # save potential trees
        # transform them into numpy arrays
        obj["trees"] = []
        for i in self.tree:
            obj["trees"].append(i.detach().cpu().numpy()) 
                    
        obj["ELBO"] = {}
        for tree in self.elbo:
            obj["ELBO"][tree] = self.elbo[tree].detach().cpu().numpy()
        
        obj["clonal_prob"] = {}
        for tree in self.clone_probs:
            obj["clonal_prob"][tree] = self.clone_probs[tree].detach().cpu().numpy()
            
        obj["potential_trees"] = []
        for tree in self.potential_trees:
            obj["potential_trees"].append(tree.detach().cpu().numpy())
        
        obj["post_params"] = self.params_post
        obj["tree_indices"] = self.tree_indices
        obj["mutations_tree"] = [self.names[i] for i in self.muts]
        obj["M"] = self.M.detach().cpu().numpy()
        obj["N"] = self.N.detach().cpu().numpy()
        obj["cell_barcode"] = self.barcodes
        obj["umapx"] = self.umapx
        obj['umapy'] = self.umapy
        obj["mutations_matrix"] = self.names
        obj["children"] = self.children
        obj["parents"] = self.parents
        
        # save posterior predictive distributions
        obj["post_predictive"] = {}
        for tree in self.post_predictive:
            obj["post_predictive"][tree] = {}
            for key in self.post_predictive[tree]:
                if key.startswith("cf") or key.startswith("obs") or key.startswith("exome"):
                    continue
                else:
                    obj["post_predictive"][tree][key] = self.post_predictive[tree][key].detach().cpu().numpy()
        
        with open(file, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
            
        print("Tree object saved as pickle!")
        
    # function to run the model on all possible trees
    def run_all_trees(self, out_dir, num_iter = 300):
        
        print("Running model for all trees")
        print("Start time: {}".format(datetime.now()))
        
        start = timer()
        
        self.sel_mutation(initial = True)

        for i in(range(len(self.names)-2)):

            self.tree = self.potential_trees
            self.tree_indices = [i for i in range(len(self.tree))]
            self.add_mutation()
            
        # run model on all trees
        self.select_tree(num_iter = num_iter, init = 50)
        
        # compute clonal assignments for the selected trees
        for i in range(len(self.tree_indices)):
            self.clonal_assignment(tree = i)
        
        # export object as pickle
        self.export_pickle(out_dir + "/" + self.name + "_out.pickle")

        # export tree class object as pickle
        with open(out_dir + "/" + self.name + "_tree.pickle", "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        end = timer()
        tim = timedelta(seconds=end-start)
        h, r = divmod(tim.seconds, 3600)
        m, s = divmod(tim.seconds, 60)
        print("Total run time {}h {}m {}s".format(h, m, s))
        
     # function to get posterior predictive distributions of latent variables (more correct than taking the parameters from the best iteration)    
    def get_post_predictive(self, tree, num_samples):
        
        pyro.clear_param_store()
        
        optimiser = pyro.optim.AdagradRMSProp({})
        # loss function
        loss_func = pyro.infer.TraceEnum_ELBO(max_plate_nesting = 1, num_particles = 5)

        # set up inference algorithm
        svi = SVI(self.model, self.guide, optimiser, loss=loss_func)
        
        # subsetting allele frequencies is different if there are several populations or not
        if(self.class_af):
            
            af_alpha = self.af_alpha[:,self.muts]
            af_beta = self.af_beta[:,self.muts]
            
        else:
            af_alpha = self.af_alpha[self.muts]
            af_beta = self.af_beta[self.muts]
        
        # iterate until lowest ELBO was obtained
        nit = self.best_iter[tree]
        
        # re-run the model for the selected tree
        for i in range(nit):

            svi.step(M = self.M[:, self.muts],
                     N = self.N[:, self.muts],
                     mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                     A = self.potential_trees[tree],
                     h_alpha = self.h_alpha,
                     h_beta = self.h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
                     bulk_af = self.bulk_af,
                     class_af = self.class_af,
                     class_names = self.class_names, 
                     class_assign = self.class_assign,
                     cnv_celltype = self.cnv_celltype,
                     celltype = self.celltype,
                     celltype_names = self.celltype_names,
                     cnv_ct_mean = self.cnv_ct_mean,
                     cnv_ct_sd = self.cnv_ct_sd)  
        
        # set Predictive class 
        pred = Predictive(self.model, guide=self.guide, num_samples=num_samples)
        
        # get samples from posterior distribution
        post = pred(M = self.M[:, self.muts],
                     N = self.N[:, self.muts],
                     mut_type = [self.mut_type.tolist()[i] for i in self.muts],
                     A = self.potential_trees[tree],
                     h_alpha = self.h_alpha,
                     h_beta = self.h_beta,
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
                     bulk_af = self.bulk_af,                    
                     class_af = self.class_af,
                     class_names = self.class_names, 
                     class_assign = self.class_assign,
                     cnv_celltype = self.cnv_celltype,
                     celltype = self.celltype,
                     celltype_names = self.celltype_names,
                     cnv_ct_mean = self.cnv_ct_mean,
                     cnv_ct_sd = self.cnv_ct_sd)
        
        # posterior predictive distributions for model evaluation
        self.post_predictive[tree] = post       
        
    # function to run the model (ideally for non-iteractive running of the mode)
    def infer_hierarchy(self, num_iter, init, out_dir, num_particles = 5, print_elbo = False):

        print("Inferring clonal hierarchies...")

        # run inference model
        for i in range(len(self.names)-1):

            if i == 0:

                # select intital two mutations based on coverage
                self.sel_mutation(initial = True)
                self.select_tree(num_iter, init)

            else:

                self.add_mutation()
                self.select_tree(num_iter, init)

        # compute clonal probabilities for selected trees
        for i in range(len(self.tree_indices)):
            self.clonal_assignment(tree = i)

        # export object as pickle
        self.export_pickle(out_dir + "/" + self.name + "_out.pickle")

        # export tree class object as pickle
        with open(out_dir + "/" + self.name + "_tree.pickle", "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

#     # function to compute posterior distribution of clonal assignment
#     def clonal_assignment_new_data(self, tree, num_iter = 0):
        
#         # run svi for the tree of interest        
#         optimiser = pyro.optim.AdagradRMSProp({})
                  
#         # compute clonal assignments
#         loss_func = pyro.infer.TraceEnum_ELBO(max_plate_nesting = 1, num_particles = 5)

#         # Infer parameters for the tree of interest
#         svi = SVI(self.model, self.guide, optimiser, loss=loss_func)
        
#         # clear previous parameters
#         pyro.clear_param_store()
        
#         if num_iter == 0:
#             num_iter = self.num_iter

#         # iterate until lowest ELBO was obtained
#         nit = self.best_iter[tree]
            
#         # run the model on the selected tree
#         for i in range(nit):
            
#             svi.step(M = self.M[:, self.muts],
#                      N = self.N[:, self.muts],
#                      mut_type = [self.mut_type.tolist()[i] for i in self.muts],
#                      A = self.tree[tree],
#                      h_alpha = self.h_alpha,
#                      h_beta = self.h_beta,
#                      af_alpha = self.af_alpha[:,self.muts],
#                      af_beta = self.af_beta[:,self.muts],
#                      r_cnv = self.r_cnv[self.muts],
#                      names = [self.names[i] for i in self.muts],
#                      class_af = self.class_af,
#                      class_names = self.class_names, 
#                      class_assign = self.class_assign,
#                      cnv_celltype = self.cnv_celltype,
#                      celltype = self.celltype,
#                      celltype_names = self.celltype_names,
#                      cnv_ct_mean = self.cnv_ct_mean,
#                      cnv_ct_sd = self.cnv_ct_sd)
        
#         print("Computing clonal assignment probabilities for tree {}".format(tree))

#         # set up inference algorithm
#         svi = SVI(self.model, self.guide_clones, optimiser, loss=loss_func)
        
#         start = timer()
        
#         for i in range(num_iter):
            
#             svi.step(M = self.M_new[:, self.muts],
#                      N = self.N_new[:, self.muts],
#                      mut_type = [self.mut_type.tolist()[i] for i in self.muts],
#                      A = self.tree[tree],
#                      h_alpha = self.h_alpha,
#                      h_beta = self.h_beta,
#                      af_alpha = self.af_alpha[:,self.muts],
#                      af_beta = self.af_beta[:,self.muts],
#                      r_cnv = self.r_cnv[self.muts],
#                      names = [self.names[i] for i in self.muts],
#                      class_af = self.class_af,
#                      class_names = self.class_names_new, 
#                      class_assign = self.class_assign_new,
#                      cnv_celltype = self.cnv_celltype,
#                      celltype = self.celltype_new,
#                      celltype_names = self.celltype_names,
#                      cnv_ct_mean = self.cnv_ct_mean,
#                      cnv_ct_sd = self.cnv_ct_sd)
    
#         # print run time
#         tree_end = timer()
#         tree_time = timedelta(seconds=tree_end-start)
#         m, s = divmod(tree_time.seconds, 60)
        
#         print("Clonal assignment probabilities computed in {}m and {}s".format(m,s))
        
#         # save clonal probabilites as attribute
#         if len(self.clone_probs) == 0:
            
#             self.clone_probs = {}
            
#         self.clone_probs[tree+1] = pyro.param("clonal_probs")
        
                    
                
                
                
            
            
            
            
