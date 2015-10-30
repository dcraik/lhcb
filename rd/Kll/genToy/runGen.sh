#!/bin/bash

if [ $# -lt 2 ]
then
    echo "Usage: $0 <nToys> <firstToy>"
fi

nToys=$1
firstToy=$2

for i in $(seq 0 $(( nToys - 1 )) )
do
    for q in {0..18}
    do
        ./gen $(( $firstToy+$i )) $q
    done
done

