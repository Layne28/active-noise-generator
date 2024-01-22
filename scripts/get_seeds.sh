#!/bin/bash

nx=$1
ny=$2

taus=(1.000000 10.000000 100.000000 1000.000000 10000.000000)
lambdas=(1.000000 3.000000 10.000000 30.000000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
	filename="seeds/seeds_nx=${nx}_ny=${ny}_tau=${tau}_lambda=${lambda}.txt"
	touch $filename
	for i in {1..5}
	do
            printf "$RANDOM\n" >> $filename
	done
    done
done
