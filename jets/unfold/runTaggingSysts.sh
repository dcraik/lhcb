#!/bin/bash

name=$1
commonOpts="--data 2016 --weight-sim 1 --dfit-skip-sumw2-fits 1 --dfit-allowMassWidthScale --dfit-allowPromptMeanShift --dfit-allowPromptWidthScale --dfit-enhanced-fits --dfit-splitBkgMassShape --dfit-splitBkgIPShape --sv-binnedtemplates --sv-ptbintemplates --sv-mistag-shape-from-back --unfold-mode bayes --unfold-reg 2 --dfit-simpleEffs --pt-corr-file" # --dfit-jetPt-binned-effs"
twodeff="--dfit-effFileD0 ../efficiencies/D0Effs_45_16x8bins_up210222.root --dfit-effFileDp ../efficiencies/DpEffs_45_16x8bins_up210220.root"

#nominal
./getTaggingEffs --dir dijets_${name} $commonOpts --dfit-use-ptfrac-bins --force-stages 3 4 5 2>&1 | tee dijets_${name}.log
#./getTaggingEffs --dir dijets_${name} $commonOpts --dfit-use-ptfrac-bins --force-stages 5 2>&1 | tee dijets_${name}.log
##
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
#cp -r dijets_${name} dijets_${name}_DSyst_combMass
#cp -r dijets_${name} dijets_${name}_DSyst_combIP
cp -r dijets_${name} dijets_${name}_DSyst_disp
cp -r dijets_${name} dijets_${name}_DSyst_disp_mean
cp -r dijets_${name} dijets_${name}_DSyst_disp_width
cp -r dijets_${name} dijets_${name}_DSyst_disp_mean_down
cp -r dijets_${name} dijets_${name}_DSyst_disp_width_down
cp -r dijets_${name} dijets_${name}_DSyst_disp_mean_up
cp -r dijets_${name} dijets_${name}_DSyst_disp_width_up
cp -r dijets_${name} dijets_${name}_DSyst_prompt
cp -r dijets_${name} dijets_${name}_DSyst_mass
cp -r dijets_${name} dijets_${name}_DSyst_noptbins
#
#####Evt-by-evt weights
./getTaggingEffs --dir dijets_${name}_evtByEvtWeights --data 2016 --weight-sim 1 --dfit-skip-sumw2-fits 0 --dfit-allowMassWidthScale --dfit-allowPromptMeanShift --dfit-allowPromptWidthScale --dfit-enhanced-fits --dfit-splitBkgMassShape --dfit-splitBkgIPShape --sv-binnedtemplates --sv-ptbintemplates --sv-mistag-shape-from-back --unfold-mode bayes --unfold-reg 2 --dfit-simpleEffs --pt-corr-file --dfit-use-ptfrac-bins --use-evt-by-evt-weighting --force-stages 3 5 2>&1 | tee dijets_${name}_evtByEvtWeights.log
#
#####SV fit systematics
./getTaggingEffs --dir dijets_${name}_svSyst_bkgrnd --data 2016 --weight-sim 1 --dfit-allowMassWidthScale --dfit-allowPromptMeanShift --dfit-allowPromptWidthScale --dfit-enhanced-fits --dfit-splitBkgMassShape --dfit-splitBkgIPShape --sv-binnedtemplates --sv-ptbintemplates --sv-mistag-shape-from-back 0 --unfold-mode bayes --unfold-reg 2 --dfit-simpleEffs --pt-corr-file --dfit-use-ptfrac-bins --force-stages 2 4 2>&1 | tee dijets_${name}_svSyst_sig.log
./getTaggingEffs --dir dijets_${name}_svSyst_sig $commonOpts --dfit-use-ptfrac-bins --sv-enhancement-passes 0 --force-stages 4 2>&1 | tee dijets_${name}_svSyst_sig.log
#
#####Jet energy
./getTaggingEffs --dir dijets_${name}_jetEnergyScaleDown $commonOpts --dfit-use-ptfrac-bins --jet-energy-scale 0.93 --force-stages 5 2>&1 | tee dijets_${name}_jetEnergyScaleDown.log
./getTaggingEffs --dir dijets_${name}_jetEnergyScaleUp $commonOpts --dfit-use-ptfrac-bins --jet-energy-scale 0.97 --force-stages 5 2>&1 | tee dijets_${name}_jetEnergyScaleUp.log
./getTaggingEffs --dir dijets_${name}_jetEnergySmearDown $commonOpts --dfit-use-ptfrac-bins --jet-energy-smear 0.11 --force-stages 5 2>&1 | tee dijets_${name}_jetEnergySmearDown.log
./getTaggingEffs --dir dijets_${name}_jetEnergySmearUp $commonOpts --dfit-use-ptfrac-bins --jet-energy-smear 0.17 --force-stages 5 2>&1 | tee dijets_${name}_jetEnergySmearUp.log
#
#####D0 fit systematics
./getTaggingEffs --dir dijets_${name}_DSyst_ptbins $commonOpts --force-stages 3 2>&1 | tee dijets_${name}_DSyst_ptbins.log
./getTaggingEffs --dir dijets_${name}_DSyst_comb $commonOpts --dfit-use-ptfrac-bins --dfit-comb-shape-syst --force-stages 3 2>&1 | tee dijets_${name}_DSyst_comb.log
#./getTaggingEffs --dir dijets_${name}_DSyst_combMass $commonOpts --dfit-use-ptfrac-bins --dfit-splitBkgMassShape 2>&1 | tee dijets_${name}_DSyst_combMass.log
#./getTaggingEffs --dir dijets_${name}_DSyst_combIP $commonOpts --dfit-use-ptfrac-bins --dfit-splitBkgMassShape --dfit-splitBkgIPShape 2>&1 | tee dijets_${name}_DSyst_combIP.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp $commonOpts --dfit-allowDisplcMeanShift --dfit-allowDisplcWidthScale --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_mean $commonOpts --dfit-allowDisplcMeanShift --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_mean.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_width $commonOpts --dfit-allowDisplcWidthScale --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_width.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_mean_down $commonOpts --dfit-shiftFixedDisplcMean -1 --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_mean_down.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_width_down $commonOpts --dfit-shiftFixedDisplcWidth -1 --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_width_down.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_mean_up $commonOpts --dfit-shiftFixedDisplcMean 1 --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_mean_up.log
./getTaggingEffs --dir dijets_${name}_DSyst_disp_width_up $commonOpts --dfit-shiftFixedDisplcWidth 1 --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_disp_width_up.log
./getTaggingEffs --dir dijets_${name}_DSyst_prompt $commonOpts --dfit-splitPromptMeanShift --dfit-splitPromptWidthScale --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_DSyst_prompt.log
./getTaggingEffs --dir dijets_${name}_DSyst_mass $commonOpts --dfit-splitMassWidthScale --dfit-use-ptfrac-bins --force-stages 3 2>&1 | tee dijets_${name}_mass.log
./getTaggingEffs --dir dijets_${name}_DSyst_noptbins $commonOpts --dfit-n-pt-bins 1 --force-stages 3 2>&1 | tee dijets_${name}_DSyst_noptbins.log

###3D PIDCalib
#./getTaggingEffs --dir dijets_${name}_3DPIDCalib $commonOpts --dfit-use-ptfrac-bins --dfit-effFileD0 ../efficiencies/D0Effs_45_16x8bins_3D_weight_up210614.root --dfit-effFileDp ../efficiencies/DpEffs_45_16x8bins_3D_weight_up210614.root --force-stages 3 5 2>&1 | tee dijets_${name}_3DPIDCalib.log

###2D PIDCalib
./getTaggingEffs --dir dijets_${name}_2DPIDCalib $commonOpts $twodeff --dfit-use-ptfrac-bins --force-stages 3 5 2>&1 | tee dijets_${name}_2DPIDCalib.log
#
#for syst in {ptbins,comb,disp,disp_mean,disp_width,prompt,mass,noptbins,disp_mean_down,disp_width_down,disp_mean_up,disp_width_up}
#do
#	root -b -q -l "printTagSyst.C(\"$name\",\"DSyst_$syst\")"
#done
