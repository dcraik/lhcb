#!/bin/bash

#dijet data
cat ../davinci/lfns/lfns577.log > lfns/lfns100.log
cat ../davinci/lfns/lfns583.log > lfns/lfns101.log
cat ../davinci/lfns/lfns633.log > lfns/lfns102.log
##NEW
cat ../davinci/lfns/lfns{25,26}.log > lfns/lfns116.log
cat ../davinci/lfns/lfns{27,28}.log > lfns/lfns117.log
cat ../davinci/lfns/lfns{29,30}.log > lfns/lfns118.log

#Z+jet data
cat ../davinci/lfns/lfns{31,32}.log > lfns/lfns105.log
cat ../davinci/lfns/lfns{19,20}.log > lfns/lfns106.log
cat ../davinci/lfns/lfns{21,22}.log > lfns/lfns107.log
cat ../davinci/lfns/lfns{23,24}.log > lfns/lfns108.log

#(mumu)SS+jet data
cat ../davinci/lfns/lfnsX{5,6}.log  > lfns/lfns135.log
cat ../davinci/lfns/lfnsX{7,8}.log  > lfns/lfns136.log
cat ../davinci/lfns/lfnsX{9,10}.log > lfns/lfns137.log
cat ../davinci/lfns/lfnsX{11,12}.log> lfns/lfns138.log

#di-charm MC
cat ../davinci/lfns/lfns{564,568}.log > lfns/lfns140.log
cat ../davinci/lfns/lfns{565,569}.log > lfns/lfns141.log
cat ../davinci/lfns/lfns{566,570}.log > lfns/lfns142.log
cat ../davinci/lfns/lfns{559,562}.log > lfns/lfns143.log
cat ../davinci/lfns/lfns{567,571}.log > lfns/lfns144.log
#di-beauty MC
cat ../davinci/lfns/lfns{596,600}.log > lfns/lfns150.log
cat ../davinci/lfns/lfns{597,601}.log > lfns/lfns151.log
cat ../davinci/lfns/lfns{598,602}.log > lfns/lfns152.log
cat ../davinci/lfns/lfns{560,563}.log > lfns/lfns153.log
cat ../davinci/lfns/lfns{599,603}.log > lfns/lfns154.log

#Z+jet MC
#cat ../davinci/lfns/lfns{5,7}.log > lfns/lfns164.log
#cat ../davinci/lfns/lfns{6,8}.log > lfns/lfns165.log
cat ../davinci/lfns/lfns{33,34}.log > lfns/lfns164.log
cat ../davinci/lfns/lfns{35,36}.log > lfns/lfns165.log

#J/psi MC
cat ../davinci/lfns/lfns{608,609}.log > lfns/lfns200.log
#D0 MC
cat ../davinci/lfns/lfns546.log > lfns/lfns201.log
#D+ MC
cat ../davinci/lfns/lfns547.log > lfns/lfns202.log
#Ds MC
cat ../davinci/lfns/lfns548.log > lfns/lfns203.log
#D0->K3pi MC
cat ../davinci/lfns/lfns549.log > lfns/lfns205.log

#di-jets again (for no-jet tagging)
cat ../davinci/lfns/lfns{564,568}.log > lfns/lfns240.log
cat ../davinci/lfns/lfns{565,569}.log > lfns/lfns241.log
cat ../davinci/lfns/lfns{566,570}.log > lfns/lfns242.log
cat ../davinci/lfns/lfns{559,562}.log > lfns/lfns243.log
cat ../davinci/lfns/lfns{567,571}.log > lfns/lfns244.log
#di-beauty MC
cat ../davinci/lfns/lfns{596,600}.log > lfns/lfns250.log
cat ../davinci/lfns/lfns{597,601}.log > lfns/lfns251.log
cat ../davinci/lfns/lfns{598,602}.log > lfns/lfns252.log
cat ../davinci/lfns/lfns{560,563}.log > lfns/lfns253.log
cat ../davinci/lfns/lfns{599,603}.log > lfns/lfns254.log
