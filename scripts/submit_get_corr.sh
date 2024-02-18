#!/bin/bash
#Submit active noise correlation 

#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=1500M

nx=$1
ny=$2
tau=$3
lambda=$4
dim=$5
seed=1

module load python/3.8.6

#mydir=${traj_file%/*}/
#filename=${traj_file%.*}.h5

mydir=${PROJECT}/active-noise/${dim}d/nx\=${nx}_ny\=${ny}/tau\=${tau}/lambda\=${lambda}/seed\=${seed}/
filename="${mydir}/noise_traj.h5"

echo $mydir
cd $mydir

echo "Current directory: "
pwd
echo "Computing correlation function with dim=${dim}"
echo $mydir
echo $filename

myscript="/jet/home/lfrechet/active-noise-generator/scripts/get_corr_${dim}d.py"

cp $myscript .

python "get_corr_${dim}d.py" $filename

