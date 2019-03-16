#!/bin/bash

f=$1
dir=$2
jobid=$3
subdir=$4

echo Downloading subjob $subdir
lb-run LHCbDIRAC dirac-dms-get-file LFN:$f -D$dir/$jobid/$subdir/ &> /dev/null
lb-run LHCbDIRAC dirac-dms-get-file LFN:${f%/*}/LumiTuple.root -D$dir/$jobid/$subdir/ &> /dev/null
echo Skimming subjob $subdir
lb-run ROOT ./skimTuples $jobid $subdir $dir &> /dev/null
echo Finishing subjob $subdir
rm -rf $dir/$jobid/$subdir/output.root
rm -rf $dir/$jobid/$subdir/LumiTuple.root
