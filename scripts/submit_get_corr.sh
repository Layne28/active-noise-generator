#!/bin/bash
#Submit active noise simulations

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=12:00:00
#SBATCH -N 1 
#SBATCH -n 1

module load share_modules/ANACONDA/5.3_py3

traj_file=$1
dim=$2

mydir=${traj_file%/*}/
cd $mydir
echo "Computing correlation function with dim=${dim}"

myscript="/home/laynefrechette/active-noise/active-noise-generator/scripts/get_corr_${dim}d.py"

python $myscript $traj_file

