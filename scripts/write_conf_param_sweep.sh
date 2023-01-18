#!/bin/bash

nx=32
ny=32
nz=32

for tau in 0.01 0.1 1.0 10.0
do
    for lambda in 1.0 2.0 4.0 8.0
    do
        python3 scripts/write_conf_file.py conf --nx $nx --ny $ny --nz $nz --tau $tau --Lambda $lambda --output_dir "/scratch0/laynefrechette/active-noise-results/nx=${nx}_ny=${ny}_nz=${nz}/tau=${tau}/lambda=${lambda}/"
    done
done