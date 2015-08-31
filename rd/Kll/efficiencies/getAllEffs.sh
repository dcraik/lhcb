#!/bin/bash

for w in {0..5}
do
  for b in {10,25,50}
  do
    echo $w $b
    root -b -q -l "getEffs.C(${w},${b})" > log/eff${w}_${b}.log
  done
done
