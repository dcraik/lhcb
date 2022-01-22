#!/bin/bash

mode=$1

for file in ${@:2}
do
	./skim.sh $mode $file&
done
