#!/bin/bash

for i in {16..18}
do
  for j in {1..100}
  do
    echo $i $j
    mkdir -p eff/${i}
    ./get_moments_B2Kstll ~/lhcb/rd/Kll/tuples/fromPatrick/Kmm_Q${i}_reduced.root eff/${i}/results_${j}.root finalTree_KMuMu  ../efficiencies/effParams_PIDonly_fixB_50.txt ${i} ../efficiencies/toys/${j}/Dveto_200.root ../efficiencies/toys/${j}/Psiveto_200.root ../massFit/bkgParams/${i}.dat 0 ${j}
  done
done
