#!/bin/bash

#. LbLogin.sh -c x86_64-slc6-gcc49-opt

#workdir="/tmp/dcraik"
#savedir="/eos/lhcb/user/d/dcraik/cepphi"
workdir="/data/cep-phi"
savedir="/data/cep-phi"

cleanup=1
rerun=1

for jobid in $@
do
	subdir=0
	
	for f in `cat ../davinci/lfns/lfns${jobid}.log`
	do
		joblist=($(jobs -p))
		while (( ${#joblist[*]} >= 12 ))
		do
			sleep 1
			joblist=($(jobs -p))
		done
		if [ "$f" != "NOFILE" ]
		then
			mkdir -p $workdir/$jobid/$subdir/
			if (( rerun )) || [ ! -f $workdir/$jobid/$subdir/skimmed.root ] || [ $workdir/$jobid/$subdir/skimmed.root -ot ./skim ] #&& [ ! -f $savedir/$jobid/$subdir/skimmed.root ]
			then
				./doSkim.sh $f $workdir $jobid $subdir $cleanup &
			fi
		else
			echo Skipping subjob $subdir, not finished
		fi
		subdir=$((subdir+1))
	done

	wait `jobs -p`

	#lb-run ROOT hadd -f $workdir/$jobid/skimmed.root $workdir/$jobid/*/skimmed.root
done
