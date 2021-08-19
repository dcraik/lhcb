#!/bin/bash

name=$1

#nominal
./getTaggingEffs --dir dijets_${name} --force-stages 3 4 2>&1 | tee dijets_${name}.log

#######avoid rerunning common stages
cp -r dijets_${name} dijets_${name}_evtByEvtWeights
cp -r dijets_${name} dijets_${name}_svSyst_sig
cp -r dijets_${name} dijets_${name}_svSyst_bkgrnd
cp -r dijets_${name} dijets_${name}_jetEnergyScaleDown
cp -r dijets_${name} dijets_${name}_jetEnergyScaleUp
cp -r dijets_${name} dijets_${name}_jetEnergySmearDown
cp -r dijets_${name} dijets_${name}_jetEnergySmearUp
cp -r dijets_${name} dijets_${name}_DSyst_ptbins
cp -r dijets_${name} dijets_${name}_DSyst_comb
cp -r dijets_${name} dijets_${name}_DSyst_disp_mean_down
cp -r dijets_${name} dijets_${name}_DSyst_disp_width_down
cp -r dijets_${name} dijets_${name}_DSyst_disp_mean_up
cp -r dijets_${name} dijets_${name}_DSyst_disp_width_up
cp -r dijets_${name} dijets_${name}_DSyst_prompt
cp -r dijets_${name} dijets_${name}_DSyst_mass
cp -r dijets_${name} dijets_${name}_DSyst_noptbins

#####Evt-by-evt weights
./getTaggingEffs --dir dijets_${name}_evtByEvtWeights --dfit-skip-sumw2-fits 1 --use-evt-by-evt-weighting --force-stages 3 2>&1 | tee dijets_${name}_evtByEvtWeights.log

#####SV fit systematics
./getTaggingEffs --dir dijets_${name}_svSyst_bkgrnd --sv-mistag-shape-from-back 0 --force-stages 2 4 2>&1 | tee dijets_${name}_svSyst_sig.log
./getTaggingEffs --dir dijets_${name}_svSyst_sig --sv-enhancement-passes 0 --force-stages 4 2>&1 | tee dijets_${name}_svSyst_sig.log

#####Jet energy
./getTaggingEffs --dir dijets_${name}_jetEnergyScaleDown --jet-energy-scale 0.93 2>&1 | tee dijets_${name}_jetEnergyScaleDown.log
./getTaggingEffs --dir dijets_${name}_jetEnergyScaleUp --jet-energy-scale 0.97 2>&1 | tee dijets_${name}_jetEnergyScaleUp.log
./getTaggingEffs --dir dijets_${name}_jetEnergySmearDown --jet-energy-smear 0.11 2>&1 | tee dijets_${name}_jetEnergySmearDown.log
./getTaggingEffs --dir dijets_${name}_jetEnergySmearUp --jet-energy-smear 0.17 2>&1 | tee dijets_${name}_jetEnergySmearUp.log

#####D0 fit systematics
./getTaggingEffs --dir dijets_${name}_DSyst_ptbins --dfit-use-ptfrac-bins 0 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_ptbins.log
./getTaggingEffs --dir dijets_${name}_DSyst_comb --dfit-comb-shape-syst --force-stages 3 2>&1 | tee dijets_${name}_DSyst_comb.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp --dfit-allowDisplcMeanShift --dfit-allowDisplcWidthScale --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_mean_down --dfit-shiftFixedDisplcMean -1 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_mean_down.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_width_down --dfit-shiftFixedDisplcWidth -1 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_width_down.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_mean_up --dfit-shiftFixedDisplcMean 1 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_mean_up.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_width_up --dfit-shiftFixedDisplcWidth 1 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_width_up.log
./getTaggingEffs --dir dijets_${name}_DSyst_prompt --dfit-splitPromptMeanShift --dfit-splitPromptWidthScale --force-stages 3 2>&1 | tee dijets_${name}_DSyst_prompt.log
./getTaggingEffs --dir dijets_${name}_DSyst_mass --dfit-splitMassWidthScale --force-stages 3 2>&1 | tee dijets_${name}_mass.log
./getTaggingEffs --dir dijets_${name}_DSyst_noptbins --dfit-n-pt-bins 1 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_noptbins.log

###2D PIDCalib
./getTaggingEffs --dir dijets_${name}_2DPIDCalib --dfit-effFileD0 ../efficiencies/D0Effs_45_16x8bins_up210222.root --dfit-effFileDp ../efficiencies/DpEffs_45_16x8bins_up210220.root --force-stages 3 2>&1 | tee dijets_${name}_2DPIDCalib.log

#for syst in {ptbins,comb,prompt,mass,noptbins,disp_mean_down,disp_width_down,disp_mean_up,disp_width_up}
#do
#	root -b -q -l "printTagSyst.C(\"$name\",\"DSyst_$syst\")"
#done
