#!/bin/bash

#$1 = nx, ny, nz

if [ "$1" == "" ]; then
    echo "Error: need system linear size"
    exit 0
fi

for tau in 0.01 0.1 1.0 10.0
do
    for lambda in 1.0 2.0 4.0 8.0
    do
	python scripts/plot_corr_chunks.py "/scratch0/laynefrechette/active-noise-results" plots $1 $lambda $tau 10
    done
done
