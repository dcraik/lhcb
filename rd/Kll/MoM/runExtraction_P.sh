#!/bin/bash

for i in {0..18}
do
  mkdir -p plots/fromPatrick/${i}
  ./get_moments_B2Kstll ~/Kll/tuples/fromPatrick/Kmm_Q${i}_reduced.root results_${i}_P.root finalTree_KMuMu  ../efficiencies/effParams_PIDonly_50.root ${i} ../efficiencies/Dveto_200.root ../efficiencies/Psiveto_200.root
  mv -f plots/*.pdf plots/fromPatrick/${i}/
done
