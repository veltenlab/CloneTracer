# MutaSeq v2

MutaSeq v2 is a methdology to add clonal resolution to scRNAseq datasets using nuclear and mitochondrial SNVs as well as CNVs. We have developed a Bayesian model which infers the clonal hierachy of the mutations in the sample and probabilistically assigns cells to clones. All required scripts and detailed explanation on how to implement the model can be found in:

* [Clonal inference](clonal_inference)

In MutaSeq v2, coverage of nuclear SNVs is increased through nested PCRs using mutation-specific primers. Scripts to design primers targeting nuclear SNVs of interest can be found in:

* [Primer design for nuclear SNVs library](primer_design)

## Contact

All scripts were written by [Sergi Beneyto Calabuig](https://www.crg.eu/en/group-members/sergi-beneyto-calabuig) a member of the [Velten lab](https://www.crg.eu/en/programmes-groups/velten-lab) at the Centre for Genomic Regulation (CRG). 

In case of problems or doubts about the scripts please [raise an issue](https://github.com/veltenlab/MutaSeq-v2/issues/new).
