#!/bin/bash
#SBATCH -J cellranger
#SBATCH -o /home/adufour/work/logs/cellranger_arc_108_output.out
#SBATCH -e /home/adufour/work/logs/cellranger_arc_108_error.out
#SBATCH --mem=32G
#SBATCH --export=ALL
##SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
#Purge any previous modules
module purge

#Load the application
module load bioinfo/bcl2fastq-2.20.0
module load bioinfo/cellranger-arc-2.0.0

# My command lines I want to run on the cluster
cd /home/adufour/work/fragencode/workspace/geneswitch/code/10x_reference_genome

cellranger-arc mkref --config=/home/adufour/save/scripts/cell_ranger_reference/sus_scrofa_arc_108.config