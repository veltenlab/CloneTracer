#!/bin/bash
#$ -N alltreesW
#$ -cwd
#$ -e log_files/alltrees_W.err
#$ -o log_files/alltrees_W.out
#$ -q gpu
#$ -l gpu=1

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load anaconda3/anaconda3
conda activate pyro

python all_trees_W.py


