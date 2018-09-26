#!/bin/bash

for jobid in $@
do
	subdir=0
	
	for f in `cat lfns/lfns${jobid}.log`
	do
		joblist=($(jobs -p))
		while (( ${#joblist[*]} >= 4 ))
		do
			sleep 1
			joblist=($(jobs -p))
		done
		if [ "$f" != "NOFILE" ]
		then
			mkdir -p /eos/user/d/dcraik/jets/$jobid/$subdir/
			if [ ! -f /eos/user/d/dcraik/jets/$jobid/$subdir/${f##*/} ]
			then
				echo Downloading subjob $subdir
				lb-run LHCbDIRAC dirac-dms-get-file LFN:$f -D/eos/user/d/dcraik/jets/$jobid/$subdir/ &
			#else
			#	echo Skipping subjob $subdir, file already exists
			fi
		else
			echo Skipping subjob $subdir, not finished
		fi
		subdir=$((subdir+1))
	done
done

wait `jobs -p`
