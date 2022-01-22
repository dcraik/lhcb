#!/bin/bash

mode=$1

for file in ${@:2}
do
	./skimMC.sh $mode $file&
done
