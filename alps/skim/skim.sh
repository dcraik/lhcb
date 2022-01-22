#!/bin/bash

mode=$1
file=$2

./TreeSkimmer $mode -1 $file &> skim${mode}_${file}.log
