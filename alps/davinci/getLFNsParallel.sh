#!/bin/bash

workdir=$1

for jobid in ${@:2}
do
	subdir=0
	
	for f in `cat lfns/lfns${jobid}.log`
	do
		joblist=($(jobs -p))
		while (( ${#joblist[*]} >= 8 ))
		do
			sleep 1
			joblist=($(jobs -p))
		done
		if [ "$f" != "NOFILE" ]
		then
			mkdir -p $workdir/$jobid/$subdir/
			echo Downloading subjob $subdir
			./getLFN.sh $f $workdir $jobid $subdir&
		else
			echo Skipping subjob $subdir, not finished
		fi
		subdir=$((subdir+1))
	done

	wait `jobs -p`
done
