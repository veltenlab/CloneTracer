# Inference of clonal hierarchies


## Description

Here we developed and bayesian model to infer clonal hierarchies from scRNAseq data using nuclear and mitochondrial SNVs as well as CNVs. 
The model uses stochastic variational inference to select the mutation tree with the highest evidence (lowest Evidence Lower Bound, ELBO).
Then it computes the posterior probability for each cell to belong to each clone in the inferred tree.


The model is implemented in [pyro](http://pyro.ai/) a probabilistic programming language written in Python which uses PyTorch as backend.  


## Installation and dependencies


The following python libraries are required to run the model: pyro, pytorch, pickle, pandas, numpy and seaborn. 

We have included a yml which contains all required packages to run the model. A conda environmet can be created from it to run the model:

```


```
