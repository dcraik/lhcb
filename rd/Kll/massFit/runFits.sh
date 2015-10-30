#!/bin/bash

for i in {0..20}
do
   root -b -q -l "fitKmm_allQ_P.C($i)" > log/fromPatrick/kmm_Q$i.log
done
