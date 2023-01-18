#!/bin/bash
#compress trajectory files

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH -n 1

tar -cvzf $1/noise_traj.tar $1/noise_*.txt --remove-files
