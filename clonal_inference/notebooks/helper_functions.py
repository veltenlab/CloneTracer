# file with helper functions to find the best tree(s) heuristically

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

# function to get all nodes in a branch       
def get_nodes_branch(node, tree):

    # get new mutations of node
    mut = torch.flatten(torch.nonzero(tree[node, :] - tree[node-1, :] == 1, as_tuple=False))

    # get nodes that belong to the branch
    nodes = torch.nonzero(tree[:, mut] == 1, as_tuple=False)

    return(nodes[:,0])

# add mutation to a node (last column is always the recently added mutation)
def mut_to_node(t, n, children):
    # list with created trees
    new_trees = []

    # number of mutations
    nmuts = t.shape[1]-1

    # create new node with the mutation of interest
    new_node = t[n, :].clone().detach()
    new_node = new_node.unsqueeze(dim=0)
    new_node[:,nmuts] = 1

    # add row with the new node
    new_tree = torch.cat((t, new_node), dim = 0)

    # if the node is root there are only 2 possibilities
    if n == 0:
        for i in range(2):
            new_trees.append(new_tree.clone().detach())
        new_trees[-1][1:,nmuts] = 1
    else:
        # add mutation in the current node 
        new_trees.append(t.clone().detach())
        new_trees[-1][n:, nmuts] = 1

        # add new trees to the list (the unmodified corresponds to the new mutation being a node indpendent from downstream mutations)
        m = 2
        if len(children[n]) == 0:
            m = 1
        for i in range(m):
            new_trees.append(new_tree.clone().detach())

        # add new mutation as downstream node 
        new_trees[-1][(n+1):,nmuts] = 1

        # check if node is a branching point and add the mutation upstream of each of the branching nodes
        if len(children[n]) > 1:
            for i in children[n]:
                new_trees.append(new_tree.clone().detach())
                brch_nodes = get_nodes_branch(i, t)
                new_trees[-1][brch_nodes, nmuts] = 1

    return(new_trees)

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
        # branches where mutations are places
        self.children = []
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
              class_af = False, class_names = [], class_assign = [], 
              cnv_celltype = False, celltype = [], celltype_names = [], cnv_ct_mean = [], cnv_ct_sd = []):

        types = ["CNV", "nuclSNV", "mitoSNV"]

        nmut = M.shape[1]
        ncells = M.shape[0]

        # create latent variables for parameter
        fpr,h,u,lnu = torch.zeros(3), torch.zeros(nmut), torch.zeros(3), torch.zeros(3)
        
        # hcnv and r can be a unique variable per chromosome or be celltype-specific
        if cnv_celltype:
            logith = torch.zeros(len(celltype_names))
            hcnv = torch.zeros(len(celltype_names))
            logitr = torch.zeros(len(celltype_names))
            r = torch.zeros(len(celltype_names))
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
                fpr[i] = pyro.sample("fpr_{}".format(n), dist.Gamma(torch.tensor(2.), torch.tensor(10.)))

            # non-informative prior for concetration parameter
            u[i] = pyro.sample("u_{}".format(n), dist.Normal(0, 3))

            # exponentiate the overdispersion parameter 
            lnu[i] = torch.exp(u[i])

        # create mutation-specific latent variables (h, u and r)
        for j in range(nmut):

            # get mutation name
            m = names[j]

            if mut_type[j] != 0:
                               
                # heteroplasmy as Beta distribution (high values e.g 1000 for nuclear SNVs, smaller for mitochondrial)
                h[j] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha[j], h_beta[j]))
                
            
            # CNVs have an additional parameter r which represent the proportion of reads in the affected chr in cancer cells
            else:
                            
                # celltype-specific CNV model
                if cnv_celltype:
                    for ind,ct in enumerate(celltype_names):
                        
                        # latent variable logit h using celltype-specific priors inferred from patients w/o CNV
                        logith[ind] = pyro.sample("hcnv_{}_{}".format(m, ct), dist.Normal(cnv_ct_mean[ind], cnv_ct_sd))
                    
                        # sigmoid transformation to compute hcnv (inverse of logit)
                        hcnv[ind] = torch.sigmoid(logith[ind])
                                                
                    # cnv_ratio in this case is the ratio between the fraction of reads in the affected chromosome between cancer and healthy
                    cnv_ratio = pyro.sample("cnv_ratio_{}".format(m), dist.Normal(r_cnv[j], 0.01))
               
                    # to compute the celltype-specific r hcvn is multiplied by r_ratio
                    r = torch.mul(hcnv, cnv_ratio)
                
                # general model which assumes that the ratio of reads in the affected chromosome is the same in all cell populations        
                else:
                    # heteroplasmy as Beta distribution (high values e.g 1000 for nuclear SNVs, smaller for mitochondrial)
                    hcnv[j] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha[j], h_beta[j]))

                    # r in this case comes from a Beta distribution. Alpha is multiplied by 1.5 or 0.5
                    r[j] = pyro.sample("r_{}".format(m), dist.Beta(h_alpha[j]*r_cnv[j], h_beta[j])) 

                    #the observation we have for this quantity is the ratio between r[j] (cancer) and hcnv[j] (healthy)
                    pyro.sample("cnv_state_{}".format(m), dist.Normal(torch.div(r[j],hcnv[j]).item(), 0.05), obs = r_cnv[j])
                    
        # loop over all cells
        with pyro.plate("cell_loop", ncells):

            # marginalise node attachment using parallel enumeration
            # this means that for each cell it iterates over all attachments in parallel. 
            node = pyro.sample("assignment", dist.Categorical(probs = cf))

            # iterate through the mutations
            for j in range(nmut):

                # get mutation name
                m = names[j]

                # transform data into tensor
                M_counts = M[:,j]
                N_counts = N[:,j]

                # make a tensor of booleans 
                k = bool_node(A, j, node)

                # get mutation type (0: CNV, 1: nuclear, 2: mitochondrial)
                t = int(mut_type[j])

                # nuclear and mitochondrial mutations
                if t != 0:             
                    # sample from the corresponding distribution 
                    pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                       dist.Poisson(fpr[t]), 
                                                                       dist.BetaBinomial(h[j]*lnu[t], (1-h[j])*lnu[t], total_count = M_counts+N_counts)), 
                                                                        obs = M_counts)
                # CNAs
                else:
                    if cnv_celltype:
                        
                        h_ct = torch.gather(hcnv, 0, celltype)
                        r_ct = torch.clamp(torch.gather(r, 0, celltype), 1e-10, 0.9999999999)
                                                              
                        # sample from beta binomial distribution
                        pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                           dist.BetaBinomial(h_ct*lnu[t], (1-h_ct)*lnu[t], total_count = N_counts),
                                                                           dist.BetaBinomial(r_ct*lnu[t], (1-r_ct)*lnu[t], total_count = N_counts)),
                                                                           obs = M_counts) 
     
                    
                    else:
                        # sample from beta binomial distribution
                        pyro.sample("obs_{}".format(m), dist.MaskedMixture(k[node], 
                                                                           dist.BetaBinomial(hcnv[j]*lnu[t], (1-hcnv[j])*lnu[t], total_count = N_counts),
                                                                           dist.BetaBinomial(r[j]*lnu[t], (1-r[j])*lnu[t], total_count = N_counts)),
                                                                           obs = M_counts) 

    # define guide
    def guide(self,M, N, A, mut_type, h_alpha, h_beta, af_alpha, af_beta, names, r_cnv, 
              class_af = False, class_names = [], class_assign = [],
              cnv_celltype = False, celltype = [], celltype_names = [], cnv_ct_mean = [], cnv_ct_sd = []):

        types = ["CNV", "nuclSNV", "mitoSNV"]

        nmut = M.shape[1]
        ncells = M.shape[0]
        nodes = A.shape[0]

        # create latent variables for each mutation
        fpr,h,u = torch.zeros(3), torch.zeros(nmut), torch.zeros(3)
        
        # hcnv and r can be a unique variable per chromosome or be celltype-specific
        if cnv_celltype:
            logith = torch.zeros(len(celltype_names))
            hcnv = torch.zeros(len(celltype_names))
            cnv_ratio = torch.zeros(len(celltype_names))
        else:
            hcnv = torch.zeros(nmut)
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

            # mean and sd of Normal distribution for concentration parameters
            u_mean_param = pyro.param("u_mean_{}".format(n), torch.tensor(0.))
            u_sd_param = pyro.param("u_sd_{}".format(n), torch.tensor(3.0), constraint = constraints.positive)                      

            # normal distribution to approximate concentration parameters
            u[i] = pyro.sample("u_{}".format(n), dist.Normal(u_mean_param, u_sd_param))

            # parameters of mito and nuclear FPR
            if i > 0:
                fpr_shape = pyro.param("fpr_shape_{}".format(n), torch.log(torch.tensor(2.)), constraint = constraints.positive)
                fpr_scale_param = pyro.param("fpr_scale_{}".format(n), torch.tensor(0.05), constraint = constraints.positive)
                fpr_scale = torch.div(1, fpr_scale_param)

                # gamma dist. to approximate posterior of fpr latent variable
                fpr = pyro.sample("fpr_{}".format(n), dist.Gamma(torch.exp(fpr_shape), fpr_scale))

        # loop over the different mutation types
        for j in range(nmut):

            # get mutation name
            m = names[j]
            
            if mut_type[j] != 0 or cnv_celltype == False:
                
                # mean and concentration of the beta distribution as parameters
                h_mean_param = pyro.param("h_mean_{}".format(m), torch.tensor(0.5))
                h_conc_param = pyro.param("h_conc_{}".format(m), torch.tensor(3.0))

                # transform to shape parameters (helps to avoid non-identifiabilities)
                h_alpha = torch.exp(h_conc_param) * torch.sigmoid(h_mean_param)
                h_beta = torch.exp(h_conc_param) * (1-torch.sigmoid(h_mean_param))

                # beta distribution to approximate posterior of heteroplasmy latent variable
                h[j] = pyro.sample("het_{}".format(m), dist.Beta(h_alpha, h_beta))
        
            # for CNVs parameters for logitr and logithcnv are created
            if mut_type[j] == 0:
                if cnv_celltype:
                    for ind, ct in enumerate(celltype_names):
                        
                        # mean and sd for celltype specific hcnv parameters
                        hcnv_mean_param = pyro.param("hcnv_mean_{}.{}".format(m, ct), cnv_ct_mean[ind])
                        hcnv_sd_param = pyro.param("hcnv_sd_{}.{}".format(m, ct), cnv_ct_sd, constraint = constraints.positive)
                        
                        # latent variable hcnv
                        logith[ind] = pyro.sample("hcnv_{}_{}".format(m, ct), dist.Normal(hcnv_mean_param, hcnv_sd_param))
                        
                    # mean and sd of Normal distribution of cnv_ratio latent variable
                    cnv_ratio_mean = pyro.param("CNVratio_mean_{}".format(m), r_cnv[j])
                    cnv_ratio_sd = pyro.param("CNVratio_sd_{}".format(m), torch.tensor(0.05), constraint = constraints.positive)

                    # latent variable cnv_ratio
                    cnv_ratio = pyro.sample("cnv_ratio_{}".format(m), dist.Normal(cnv_ratio_mean, cnv_ratio_sd))
                        
                                            
                else:
                    # mean and concentration of the beta distribution as parameters
                    r_mean_param = pyro.param("r_mean_{}".format(m), torch.tensor(0.5))
                    r_conc_param = pyro.param("r_conc_{}".format(m), torch.tensor(0.0))

                    # beta distribution to approximate posterior of heteroplasmy latent variable
                    r_alpha = torch.exp(r_conc_param) * torch.sigmoid(r_mean_param)
                    r_beta = torch.exp(r_conc_param) * (1-torch.sigmoid(r_mean_param))

                    # latent variable r
                    r[j] = pyro.sample("r_{}".format(m), dist.Beta(r_alpha, r_beta))  
                
    # define guide for clonal assignment
    @config_enumerate
    def guide_clones(self, M, N, A, mut_type, h_alpha, h_beta, af_alpha, af_beta, names, r_cnv, 
                     class_af = False, class_names = [], class_assign = [],
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
    def select_tree(self, num_iter, init, num_particles = 5, vect_particles=False, print_elbo = False):
        
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
                                         h_alpha = self.h_alpha[self.muts],
                                         h_beta = self.h_beta[self.muts],
                                         af_alpha = af_alpha,
                                         af_beta = af_beta,
                                         r_cnv = self.r_cnv[self.muts],
                                         names = [self.names[i] for i in self.muts],
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
                
        # add tree structure as attribute
        self.get_children()
        
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
        if print_elbo:
            self.print_elbo(num_iter = num_iter, init = init)
                   
    # add mutation to selected tree 
    def add_mutation(self):
        
        # list of new trees
        new_trees = []
        
        # iterate over selected trees 
        for ind,t in enumerate(self.tree):
            # add column for the new mutation
            t = torch.cat((t, torch.zeros(t.shape[0],1)), dim =1)
            
            # iterate through the nodes and attach new mutation in all possible sites
            for n in range(t.shape[0]):
                
                new_trees.append(mut_to_node(t, n, self.children[ind]))
                
        # order nodes by number of mutations and store in potential_trees slot
        self.potential_trees = [i[torch.argsort(torch.sum(i, dim = 1)), :] for l in new_trees for i in l]
        
        # choose mutation to add based on coverage
        mut = self.sel_mutation()
        
        print("{} added to the tree".format(self.names[mut]))
            
    # function to get the children of each node
    def get_children(self):
        
        self.children = []
        
        # iterate over the trees
        for ind,t in enumerate(self.tree):

            # set children
            self.children.append([ [] for _ in range(t.shape[0])])

            # sort nodes by number of mutations
            t = t[torch.argsort(torch.sum(t, dim = 1)), :]

            num_b = 0
            # For each branch the 1st integer indicates the last mutation of the terminal node and the 2nd item refers to the index of the terminal node
            branches = {}

            # iterate through the nodes
            for i in range(t.shape[0]):

                # for root the immediate node is attached
                if i == 0:
                    self.children[ind][i] = [1]
                    continue

                # get the node
                n = t[i,]

                # for the 1st node make a new branch and attach it
                if i == 1:
                    last_mut = int(torch.nonzero(n == 1)[0])
                    branches = {0: [(last_mut,1)]}
                    num_b += +1
                    continue

                attached = False
                #iterate through branches of the tree
                for b in branches:

                    #get terminal node of branch
                    parent = t[branches[b][-1][1],]
                    p_mut = branches[b][-1][0]

                    # check whether node is attached to the branch
                    if n[p_mut] == parent[p_mut]:

                        # get the new mutation of the current node
                        last_mut = int(torch.nonzero(n-parent == 1)[0])
                        # add node to branch
                        branches[b].append((last_mut, i)) 

                        # add current node as children of the previous node
                        self.children[ind][i-1] = [i]

                        attached = True

                        break

                b_copy = branches.copy()

                #iterate through branches and nodes to determine where does the branching node attach
                if not attached:
                    for b in b_copy:
                        
                        if not attached:

                            for index, f in enumerate(reversed(range(len(b_copy[b])))):


                                j = b_copy[b][f]

                                parent = t[j[1],]
                                p_mut = j[0]

                                if n[p_mut] == parent[p_mut]:
                                    # add node as children of stem node
                                    att_n = j[1]
                                    self.children[ind][att_n].append(i)

                                    last_mut = int(torch.nonzero(n-parent)[0])

                                    # create a new branch
                                    branches[num_b] = branches[b][0:f+1]
                                    branches[num_b].append((last_mut,i))
                                    num_b += 1
                                    attached = True
                                    break

                    # if no attachment is found then it is attached to the root
                    if not attached:
                        self.children[ind][0].append(i)
                        # create a new branch
                        last_mut = int(torch.nonzero(n-torch.zeros(n.shape[0]))[0])
                        branches[num_b] = [(last_mut,i)]
                        num_b += 1
                        attached = True
                                    
        return(self.children)
    
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


        fig, ax = plt.subplots(nrows = 8, ncols = 5, figsize = (35, 50))
        count = 0
        
        if it == 0:
            for index, p in enumerate(params):

                count += 1
                row = divmod(count, 5)[0]
                col = divmod(count, 5)[1]

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
                row = divmod(count, 5)[0]
                col = divmod(count, 5)[1]
                
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
            
            fig, ax = plt.subplots(nrows = 1, ncols = num_plots, figsize = 
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
                     h_alpha = self.h_alpha[self.muts],
                     h_beta = self.h_beta[self.muts],
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
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
                     h_alpha = self.h_alpha[self.muts],
                     h_beta = self.h_beta[self.muts],
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
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
            self.get_children()
            self.add_mutation()
            
        # run model for all trees    
        self.select_tree(num_iter = num_iter, init = 50, print_elbo = False)
        
        # compute clonal probabilities for selected trees
        for i in range(len(self.tree_indices)):
            self.clonal_assignment(tree = i)
        
        # save output as pickle
        self.export_pickle(out_dir + "/all_trees_" + self.name + ".pickle")
        
        # save tree object as pickle
        with open(out_dir + "/all_trees_" + self.name + "_tree.pickle", "wb") as f:
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
                     h_alpha = self.h_alpha[self.muts],
                     h_beta = self.h_beta[self.muts],
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
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
                     h_alpha = self.h_alpha[self.muts],
                     h_beta = self.h_beta[self.muts],
                     af_alpha = af_alpha,
                     af_beta = af_beta,
                     r_cnv = self.r_cnv[self.muts],
                     names = [self.names[i] for i in self.muts],
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


        
                    
                
                
                
            
            
            
            