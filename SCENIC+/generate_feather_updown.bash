#!/bin/bash
#SBATCH -J cistarget_db
#SBATCH -c 20
#SBATCH --mem=180G
#SBATCH --constraint=broadwell
#SBATCH --output=/home/adufour/work/logs/cistarget_db.log

#initialise l environnement et les variables
version=Sus_scrofa.Sscrofa11.1.108
fasta=/home/adufour/work/cistargetdb/file_feather_generation/fasta_gene_ensembl_updown_$version.fasta
motif_dir=~/work/cistargetdb/singletons_v2
motif2tf=/home/adufour/work/cistargetdb/motif2tf_orthologuous.tbl

#recuperation des motifs
ls $motif_dir | sed -e 's/.cb//' > ~/work/cistargetdb/motif_list.txt

source /home/adufour/.bashrc
conda activate create_cistarget_databases
create_cistarget_databases_dir="/home/adufour/work/programms/create_cisTarget_databases"

/home/adufour/work/programms/create_cisTarget_databases/create_cistarget_motif_databases.py --fasta $fasta \
--motifs_dir $motif_dir \
--motifs ~/work/cistargetdb/motif_list.txt \
--genes "#[0-9]+$" \
--threads 20 \
--output ~/work/cistargetdb/feather_v3_updown/Sus_scrofa