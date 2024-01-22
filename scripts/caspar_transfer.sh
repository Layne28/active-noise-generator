#!/bin/bash
#Transfer data to caspar

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=72:00:00
#SBATCH -N 1 
#SBATCH -n 1

module load share_modules/ANACONDA/5.3_py3

from_folder=$1
to_folder=$2

echo "Transferring data from $from_folder to $to_folder"

rsync -ravz $from_folder laynefrechette@caspar.brandeis.edu:/current/laynefrechette/$to_folder/
