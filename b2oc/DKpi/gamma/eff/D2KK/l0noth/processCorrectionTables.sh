#!/bin/bash

for file in `ls tables2012L0HadronTOS_S20realET/tables_S20`
do
   python -b cleanFile.py $file
   root -b -q -l "makeTableTree.C(\"$file\")"
   root -b -q -l "makeTableHist.C(\"$file\")"
#   echo $file
#   sed 's/[^0-9\.\-]\+/\t/g' < tables_S20/$file | cut -f2,4,5,6,7 > $file.dat
done
