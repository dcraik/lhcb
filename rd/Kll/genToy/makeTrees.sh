#!/bin/bash

if [ $# -lt 3 ]
then
    echo "Usage: $0 <nToys> <firstToy> <dir>"
fi

nToys=$1
firstToy=$2
dir=$3

for i in $(seq 0 $(( nToys - 1 )) )
do
    for q in {0..18}
    do
        root -b -q -l "makeTree.C (\"$dir\",$(( $firstToy+$i )),$q)"
    done
done

