#!/bin/bash
#SBATCH -p workq
#SBATCH --mem=50G
#SBATCH -J ARC
#SBATCH -N 1
#SBATCH -o /home/adufour/work/logs/ARC_output.out
#SBATCH -e /home/adufour/work/logs/ARC_error.out
#Purge any previous modules
module purge

#Load the application
module load bioinfo/bcl2fastq-2.20.0
module load bioinfo/cellranger-arc-2.0.0

#mkdir ~/work/run_cellranger_count
cd /home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl

cellranger-arc count --id=arc_lw7_108 \
                     --reference=/home/adufour/work/10x_reference_genome/sus_scrofa_11_ensembl_arc_108 \
                     --libraries=/home/adufour/work/fragencode/workspace/plus4pigs/data/reads/LW7.1_libraries_sheets.csv \
                     --jobmode=slurm \
                     --mempercore=8 \
                     --maxjobs=100

cellranger-arc count --id=arc_lw9_108 \
                     --reference=/home/adufour/work/10x_reference_genome/sus_scrofa_11_ensembl_arc_108 \
                     --libraries=/home/adufour/work/fragencode/workspace/plus4pigs/data/reads/LW9.1_libraries_sheets.csv \
                     --jobmode=slurm \
                     --mempercore=8 \
                     --maxjobs=100