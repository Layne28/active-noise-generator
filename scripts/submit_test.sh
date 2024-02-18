#!/bin/bash
#Submit active noise correlation 

#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

module load python/3.8.6

cd ${PROJECT}/active-noise/2d/nx\=200_ny\=200/tau\=1000.000000/lambda\=10.000000/seed\=1/
pwd
touch "test.txt"
