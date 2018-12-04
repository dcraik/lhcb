#!/bin/bash

~/UraniaDev_v7r0/run python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -b calibBinnings.py -s pions2 "Turbo16" "MagDown" "Pi" "[Brunel_MC15TuneV1_ProbNNpi > 0.1]" "Brunel_P" "Brunel_PT" ""
~/UraniaDev_v7r0/run python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -b calibBinnings.py -s pions2 "Turbo16" "MagUp" "Pi" "[Brunel_MC15TuneV1_ProbNNpi > 0.1]" "Brunel_P" "Brunel_PT" ""
~/UraniaDev_v7r0/run python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -b calibBinnings.py -s kaons4 "Turbo16" "MagDown" "K" "[Brunel_MC15TuneV1_ProbNNK > 0.2]" "Brunel_P" "Brunel_PT" ""
~/UraniaDev_v7r0/run python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -b calibBinnings.py -s kaons4 "Turbo16" "MagUp" "K" "[Brunel_MC15TuneV1_ProbNNK > 0.2]" "Brunel_P" "Brunel_PT" ""
