#!/bin/bash

#SetupProject LHCbDirac v8r2p55
#~/UraniaDev_v7r0/run bash --norc

#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s pions  "Turbo16" "MagDown" "Pi" "[Brunel_MC15TuneV1_ProbNNpi > 0.2]"
#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s pions  "Turbo16" "MagUp" "Pi" "[Brunel_MC15TuneV1_ProbNNpi > 0.2]"

#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s kaons3 "Turbo16" "MagDown" "K" "[Brunel_MC15TuneV1_ProbNNK > 0.3]"
#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s kaons3 "Turbo16" "MagUp" "K" "[Brunel_MC15TuneV1_ProbNNK > 0.3]"

python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s protons6 "Turbo16" "MagDown" "P" "[Brunel_MC15TuneV1_ProbNNp > 0.3]"
python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s protons6 "Turbo16" "MagUp" "P" "[Brunel_MC15TuneV1_ProbNNp > 0.3]"

#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s misIDAsMu -o "muons/" "Turbo16" "MagDown" "Pi" "[IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5]"
#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s misIDAsMu -o "muons/" "Turbo16" "MagDown" "K" "[IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5]"
#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s misIDAsMu -o "muons/" "Turbo16" "MagDown" "P" "[IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5]"
#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s muons5 -o "muons/" "Turbo16" "MagDown" "Mu" "[IsMuon == 1 && Brunel_MC15TuneV1_ProbNNmu > 0.5]"
#python ~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py -X "Brunel_P" -Y "Brunel_PT" -Z "" -b DKpiBinning.py -s misIDAsMu -o "muons/" "Turbo16" "MagDown" "e" "[Brunel_MC15TuneV1_ProbNNmu > 0.5]"
