#!/bin/bash

input=20210626 #...08 #...0304 #...0222 #...17 #20210131_rerun
output=20210627 #...09 #...0305 #...0227 #...17 #20210131_rerun2

#baseline
./fitZJetSVs --data 201X --dir zjets_${output} --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-correct --sv-mistag-shape-from-back 1 --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 --tag-eff-file-is-root 0 --tag-eff-file jetTagComb.txt 2>&1 | tee zjets_${output}.log
#correct mistag
#./fitZJetSVs --data 201X --dir zjets_${output}_correctMistag_noTemplates --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --sv-light-yield-correct --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_correctMistag_noTemplates.log
##float mistag
./fitZJetSVs --data 201X --dir zjets_${output}_floatMisTag --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --sv-light-yield-float 1 --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 --tag-eff-file-is-root 0 --tag-eff-file jetTagComb.txt 2>&1 | tee zjets_${output}_floatMisTag.log
##jetEnergyScaleDown
./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergyScaleDown --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergyScaleDown --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-scale 0.93 2>&1 | tee zjets_${output}_jetEnergyScaleDown.log
##jetEnergyScaleUp
./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergyScaleUp --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergyScaleUp --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-scale 0.97 2>&1 | tee zjets_${output}_jetEnergyScaleUp.log
##jetEnergySmearDown
./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergySmearDown --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergySmearDown --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-smear 0.11 2>&1 | tee zjets_${output}_jetEnergySmearDown.log
##jetEnergySmearUp
./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergySmearUp --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergySmearUp --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-smear 0.17 2>&1 | tee zjets_${output}_jetEnergySmearUp.log
##svSyst_bkgrnd
./fitZJetSVs --data 201X --dir zjets_${output}_svSyst_bkgrnd --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 0 --input-dir dijets_${input}_svSyst_bkgrnd --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_svSyst_bkgrnd.log
##svSyst_sig
./fitZJetSVs --data 201X --dir zjets_${output}_svSyst_sig --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_svSyst_sig --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_svSyst_sig.log
##unfoldYBins
./fitZJetSVs --data 201X --dir zjets_${output}_unfoldYBins --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 --tag-eff-file-is-root 0 --tag-eff-file jetTagComb.txt --unfold-rapiditybins 2>&1 | tee zjets_${output}_unfoldYBins.log

##D fits and weights - these only differ from the baseline in the efficiency but this is the easiest way to produce the required files
#for syst in {DSyst_comb,DSyst_prompt,DSyst_disp_width_up,DSyst_disp_width_down,DSyst_disp_mean_up,DSyst_disp_mean_down,DSyst_mass,DSyst_noptbins,evtByEvtWeights}
#do
#    cp -r zjets_${output} zjets_${output}_${syst}
#    ./fitZJetSVs --data 201X --dir zjets_${output}_${syst} --sv-binnedtemplates --unfold-reg 2 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-correct --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_${syst} --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_${syst}.log
#done

##baseline
#./fitZJetSVs --data 201X --dir zjets_${output} --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}.log
##correct mistag
##./fitZJetSVs --data 201X --dir zjets_${output}_correctMistag_noTemplates --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --sv-light-yield-correct --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_correctMistag_noTemplates.log
###float mistag
#./fitZJetSVs --data 201X --dir zjets_${output}_floatMisTag --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --sv-light-yield-float 1 --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_floatMisTag.log
###jetEnergyScaleDown
#./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergyScaleDown --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergyScaleDown --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-scale 0.93 2>&1 | tee zjets_${output}_jetEnergyScaleDown.log
###jetEnergyScaleUp
#./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergyScaleUp --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergyScaleUp --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-scale 0.97 2>&1 | tee zjets_${output}_jetEnergyScaleUp.log
###jetEnergySmearDown
#./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergySmearDown --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergySmearDown --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-smear 0.11 2>&1 | tee zjets_${output}_jetEnergySmearDown.log
###jetEnergySmearUp
#./fitZJetSVs --data 201X --dir zjets_${output}_jetEnergySmearUp --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_jetEnergySmearUp --sv-minmcor 500 --sv-nmcorbins 38 --jet-energy-smear 0.17 2>&1 | tee zjets_${output}_jetEnergySmearUp.log
###svSyst_bkgrnd
#./fitZJetSVs --data 201X --dir zjets_${output}_svSyst_bkgrnd --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 0 --input-dir dijets_${input}_svSyst_bkgrnd --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_svSyst_bkgrnd.log
###svSyst_sig
#./fitZJetSVs --data 201X --dir zjets_${output}_svSyst_sig --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-mistag-shape-from-back 1 --input-dir dijets_${input}_svSyst_sig --sv-minmcor 500 --sv-nmcorbins 38 2>&1 | tee zjets_${output}_svSyst_sig.log
###unfoldYBins
#./fitZJetSVs --data 201X --dir zjets_${output}_unfoldYBins --sv-binnedtemplates --unfold-reg 1 --pt-bins 15e3 20e3 30e3 50e3 100e3 --y-bins 2. 2.75 3.5 4.5 --pt-corr-file --sv-light-yield-from-file zjets_20210305_correctMistag/mistagYieldConstraints.root --sv-mistag-shape-from-back 1 --input-dir dijets_${input} --sv-minmcor 500 --sv-nmcorbins 38 --unfold-rapiditybins 2>&1 | tee zjets_${output}_unfoldYBins.log
