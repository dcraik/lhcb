#!/bin/bash

#. LbLogin.sh -c x86_64-slc6-gcc49-opt

workdir=$1
savedir=$1

cleanup=0
rerun=0

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
			if (( rerun )) || [ ! -f $workdir/$jobid/$subdir/skimmed.root ] || [ $workdir/$jobid/$subdir/skimmed.root -ot ./skimTuples ] #&& [ ! -f $savedir/$jobid/$subdir/skimmed.root ]
			then
				./doSkim.sh $f $workdir $jobid $subdir $cleanup&
			fi
		else
			echo Skipping subjob $subdir, not finished
		fi
		subdir=$((subdir+1))
	done

	wait `jobs -p`

	#lb-run ROOT hadd -f $workdir/$jobid/skimmed.root $workdir/$jobid/*/skimmed.root
done
