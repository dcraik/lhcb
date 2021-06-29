#!/bin/bash

input=20210628_test #...6 #...08 #...0304 #...0222 #...17 #20210131_rerun
output=20210628_test #...7 #...09 #...0305 #...0227 #...17 #20210131_rerun2

#baseline
./fitZJetSVs --dir zjets_${output} --input-dir dijets_${input} --tag-eff-file-is-root 0 --tag-eff-file jetTagComb.txt 2>&1 | tee zjets_${output}.log
#correct mistag
#./fitZJetSVs --dir zjets_${output}_correctMistag_noTemplates --sv-light-yield-correct --input-dir dijets_${input} 2>&1 | tee zjets_${output}_correctMistag_noTemplates.log
##float mistag
./fitZJetSVs --dir zjets_${output}_floatMisTag --sv-light-yield-float 1 --input-dir dijets_${input} --tag-eff-file-is-root 0 --tag-eff-file jetTagComb.txt 2>&1 | tee zjets_${output}_floatMisTag.log
##jetEnergyScaleDown
./fitZJetSVs --dir zjets_${output}_jetEnergyScaleDown --input-dir dijets_${input}_jetEnergyScaleDown --jet-energy-scale 0.93 2>&1 | tee zjets_${output}_jetEnergyScaleDown.log
##jetEnergyScaleUp
./fitZJetSVs --dir zjets_${output}_jetEnergyScaleUp --input-dir dijets_${input}_jetEnergyScaleUp --jet-energy-scale 0.97 2>&1 | tee zjets_${output}_jetEnergyScaleUp.log
##jetEnergySmearDown
./fitZJetSVs --dir zjets_${output}_jetEnergySmearDown --input-dir dijets_${input}_jetEnergySmearDown --jet-energy-smear 0.11 2>&1 | tee zjets_${output}_jetEnergySmearDown.log
##jetEnergySmearUp
./fitZJetSVs --dir zjets_${output}_jetEnergySmearUp --input-dir dijets_${input}_jetEnergySmearUp --jet-energy-smear 0.17 2>&1 | tee zjets_${output}_jetEnergySmearUp.log
##svSyst_bkgrnd
./fitZJetSVs --dir zjets_${output}_svSyst_bkgrnd --sv-mistag-shape-from-back 0 --input-dir dijets_${input}_svSyst_bkgrnd 2>&1 | tee zjets_${output}_svSyst_bkgrnd.log
##svSyst_sig
./fitZJetSVs --dir zjets_${output}_svSyst_sig --input-dir dijets_${input}_svSyst_sig 2>&1 | tee zjets_${output}_svSyst_sig.log
##unfoldYBins
./fitZJetSVs --dir zjets_${output}_unfoldYBins --input-dir dijets_${input} --tag-eff-file-is-root 0 --tag-eff-file jetTagComb.txt --unfold-rapiditybins 2>&1 | tee zjets_${output}_unfoldYBins.log
