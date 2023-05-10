# CloneTracer

CloneTracer is a methdology to add clonal resolution to scRNAseq datasets using nuclear and mitochondrial SNVs as well as CNVs. Coverage in sample-specific nuclear SNVs as well as the mitochondrial genome is increased by generating separate targeted libraries from the whole transcriptome cDNA library (Optimized 10x libraries). 

<p align="center">
<img src="method_cartoon.png" width="400" height="350">
</p>

For details of the method see our manuscript:

[Beneyto-Calabuig, Merbach et al,. Clonally resolved single-cell multi-omics identifies routes of cellular differentiation in acute myeloid leukemia, Cell Stem Cell 2023](external.ink?to=https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(23)00119-4)
  
To analyse the generated data, we have developed CloneTracer, a Bayesian framework which infers the mutational hierarchy in a particular sample and probabilistically assigns single cells to clones. All required scripts and detailed explanation on how to implement the model can be found in:

* [Clonal inference](clonal_inference)

In order to increase the coverage on nuclear SNVs a series of nested PCR are carried out on the whole transcriptome library. This requires mutation-specific primers. All necessary scripts to design such primers are available in:

* [Primer design for nuclear SNVs library](primer_design)

Workflows to pre-process the mitochondrial and nuclear SNV Optimized 10x libraries can be found in:

* [Processing pipeline nuclear SNV library](library_processing/nuclear-snv)
* [Processing pipeline mitochondrial library](library_processing/mitochondria)

## Contact

All scripts were written by [Sergi Beneyto Calabuig](https://www.crg.eu/en/group-members/sergi-beneyto-calabuig) a member of the [Velten lab](https://www.crg.eu/en/programmes-groups/velten-lab) at the Centre for Genomic Regulation (CRG) in Barcelona. 

In case of problems or doubts about the scripts please [raise an issue](https://github.com/veltenlab/CloneTracer/issues/new).
