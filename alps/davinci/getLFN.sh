#!/bin/bash

f=$1
dir=$2
jobid=$3
subdir=$4

if [ ! -f $dir/$jobid/$subdir/output.root ]
then
	lb-run LHCbDIRAC dirac-dms-get-file LFN:$f -D$dir/$jobid/$subdir/ &> /dev/null
fi
if [ ! -f $dir/$jobid/$subdir/LumiTuple.root ]
then
	lb-run LHCbDIRAC dirac-dms-get-file LFN:${f%/*}/LumiTuple.root -D$dir/$jobid/$subdir/ &> /dev/null
fi
