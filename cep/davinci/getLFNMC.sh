#!/bin/bash

f=$1
dir=$2
jobid=$3
subdir=$4

if [ ! -f $dir/$jobid/$subdir/Tuples.root ]
then
	echo Downloading subjob $subdir
	lb-run LHCbDIRAC dirac-dms-get-file LFN:$f -D$dir/$jobid/$subdir/ &> /dev/null
fi
