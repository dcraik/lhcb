#!/bin/bash

for i in {0..20}
do
   root -b -q -l "fitKmm_allQ_P_floatRatios.C($i)" > log/fromPatrick/floatRatios/kmm_Q$i.log
done
