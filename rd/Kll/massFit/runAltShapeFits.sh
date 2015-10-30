#!/bin/bash

for i in {0..20}
do
   root -b -q -l "fitKmm_allQ_P_altShapes.C($i)" > log/fromPatrick/altShapes/kmm_Q$i.log
done
