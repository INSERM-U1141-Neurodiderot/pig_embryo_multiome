#!/bin/bash
#SBATCH -J scenic_plus
#SBATCH -c 30
#SBATCH --mem=700G
#SBATCH --output=/home/adufour/work/logs/scenic_plus_embryos.log

module purge

source /home/adufour/.bashrc
source activate scenicplus

python /home/adufour/save/scripts/scenicplus/omicscenic.py