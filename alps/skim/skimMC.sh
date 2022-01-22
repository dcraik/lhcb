#!/bin/bash

mode=$1
file=$2

#./SimSkimmer $mode -1 $file &> skim${mode}_${file}.log
./EffSkimmer $mode -1 $file &> skim${mode}_${file}.log
