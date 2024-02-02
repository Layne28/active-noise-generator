#!/bin/bash
#Submit active noise simulations

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=72:00:00
#SBATCH -N 1 
#SBATCH -n 1

module load share_modules/ANACONDA/5.3_py3

input_file=$1
nsteps=$2
seed=$3

echo "Running active noise generator with input file: $input_file"

run_dir="/home/laynefrechette/active-noise/active-noise-generator/bin/"

$run_dir/active_noise_generator $input_file $nsteps $seed

