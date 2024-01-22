#!/bin/bash

nx='200'
ny='200'

for tau in 1.000000 10.000000 100.000000 1000.000000 10000.000000
do
    for lambda in 1.000000 3.000000 10.000000 30.000000
    do
        python3 scripts/write_input_file.py input_files --nx $nx --ny $ny --tau $tau --Lambda $lambda --output_dir "/scratch0/laynefrechette/active-noise-results/2d/nx=${nx}_ny=${ny}/tau=${tau}/lambda=${lambda}/"
    done
done
