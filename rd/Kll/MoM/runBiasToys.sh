#!/bin/bash

for q in {13..16}
do
  mkdir -p toy/${q}
  for i in {0..99}
  do
    ./get_moments_B2Kstll ~/lhcb/rd/Kll/genToy/toy/toy_${i}_${q}.root toy/${q}/results_${i}.root toy  ../efficiencies/effParams_PIDonly_fixB_50.txt ${q} ../efficiencies/Dveto_200.root ../efficiencies/Psiveto_200.root "" 0
  done
done
