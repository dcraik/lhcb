#!/bin/bash

f=$1
dir=$2
jobid=$3
subdir=$4
cleanup=$5

if [ ! -f $dir/$jobid/$subdir/Tuples.root ]
then
	echo Downloading subjob $subdir
	lb-run LHCbDIRAC dirac-dms-get-file LFN:$f -D$dir/$jobid/$subdir/ &> /dev/null
fi
echo Skimming subjob $subdir
lb-run ROOT ./skim $jobid $subdir $dir &> /dev/null
echo Finishing subjob $subdir
if (( $cleanup ))
then
	rm -rf $dir/$jobid/$subdir/Tuples.root
fi
