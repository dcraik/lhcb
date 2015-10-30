#!/bin/bash

for i in {0..18}
do
  mkdir -p plots/fromPatrick/floatRatios/${i}
  ./get_moments_B2Kstll ~/lhcb/rd/Kll/tuples/fromPatrick/Kmm_Q${i}_reduced.root results_${i}_P_floatRatios.root finalTree_KMuMu  ../efficiencies/effParams_PIDonly_fixB_50.txt ${i} ../efficiencies/Dveto_200.root ../efficiencies/Psiveto_200.root ../massFit/bkgParams/floatRatios/${i}.dat
  mv -f plots/*.pdf plots/fromPatrick/floatRatios/${i}/
done
