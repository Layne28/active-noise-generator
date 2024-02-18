#!/bin/bash
#Submit active noise simulations

#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=16G

input_file=$1
nsteps=$2
seed=$3

echo "Running active noise generator with input file: $input_file"

run_dir="/jet/home/lfrechet/.local/bin/"

$run_dir/active_noise_generator $input_file $nsteps $seed

