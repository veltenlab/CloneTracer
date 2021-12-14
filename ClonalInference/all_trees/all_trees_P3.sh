#!/bin/bash
#$ -N alltreesA
#$ -cwd
#$ -e log_files/alltrees_A.err
#$ -o log_files/alltrees_A.out
#$ -q gpu_long
#$ -l gpu=1
#$ -l virtual_free=60G
#$ -l h_rt=72:00:00

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load anaconda3/anaconda3
conda activate pyro

python all_trees_A.py


