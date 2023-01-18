#!/bin/bash

#$1 = nx, ny, nz

if [ "$1" == "" ]; then
    echo "Error: need system linear size"
    exit 0
fi

rm -r "/scratch0/laynefrechette/active-noise-results/nx=${1}_ny=${1}_nz=${1}/"
