#!/bin/bash
#Submit active noise simulations

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=24:00:00
#SBATCH -N 1 
#SBATCH -n 1

module load share_modules/ANACONDA/5.3_py3

conf_file=$1

echo "Running analysis with conf file: $conf_file"

run_dir="/home/laynefrechette/active-noise/active-noise-generator/bin/"

$run_dir/analyze $conf_file

