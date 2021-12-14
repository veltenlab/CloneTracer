#!/bin/bash
#$ -N alltreesK
#$ -cwd
#$ -e log_files/alltrees_K.err
#$ -o log_files/alltrees_K.err
#$ -q gpu_long
#$ -l gpu=1
#$ -l virtual_free=60G
#$ -l h_rt=72:00:00

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load anaconda3/anaconda3
conda activate pyro

python all_trees_K.py


