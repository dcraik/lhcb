CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	skim_for_yandex	skim_for_yandex_data skim_for_yandex_data_z0 skim_for_yandex_data_light skim_for_yandex_DMC skimD0s skimD0sFromB skimCharm charmEfficiencies weightCharmEffs fit fitD0 fitD fitDs fitLc fit2d fit2dToy tagEvents skimZs fitD0withZ getTotalCharmYields getTotalCharmYieldsNew getTotalCharmYieldsMC getTotalCharmYieldsDMC getTotalCharmYieldsJMC makeNewCalibTuples makeNewDTuples fitCharmNew jpsiEff

velo.o : velo.cc velo.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex.o : skim_for_yandex.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex : skim_for_yandex.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data.o : skim_for_yandex_data.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data : skim_for_yandex_data.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data_backwards.o : skim_for_yandex_data_backwards.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data_backwards : skim_for_yandex_data_backwards.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data_new.o : skim_for_yandex_data_new.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data_new : skim_for_yandex_data_new.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data_z0.o : skim_for_yandex_data_z0.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data_z0 : skim_for_yandex_data_z0.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data_light.o : skim_for_yandex_data_light.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data_light : skim_for_yandex_data_light.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_DMC.o : skim_for_yandex_DMC.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_DMC : skim_for_yandex_DMC.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skimD0s.o : skimD0s.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimD0s : skimD0s.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skimD0sFromB.o : skimD0sFromB.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimD0sFromB : skimD0sFromB.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skimCharm.o : skimCharm.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimCharm : skimCharm.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

charmEfficiencies.o : charmEfficiencies.C charmEfficiencies.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

charmEfficiencies : charmEfficiencies.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

charmEfficienciesOld.o : charmEfficienciesOld.C charmEfficienciesOld.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

charmEfficienciesOld : charmEfficienciesOld.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

weightCharmEffs.o : weightCharmEffs.C weightCharmEffs.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

weightCharmEffs : weightCharmEffs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fit.o : fit.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fit : fit.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitD0.o : fitD0.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitD0 : fitD0.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitD.o : fitD.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitD : fitD.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitDs.o : fitDs.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitDs : fitDs.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitLc.o : fitLc.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitLc : fitLc.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fit2d.o : fit2d.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fit2d : fit2d.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fit2dnew.o : fit2dnew.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fit2dnew : fit2dnew.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fit2dToy.o : fit2dToy.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fit2dToy : fit2dToy.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fit1dnew.o : fit1dnew.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fit1dnew : fit1dnew.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

RooPromptShape.o : RooPromptShape.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

tagEvents.o : tagEvents.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

tagEvents : tagEvents.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

tagEventsNew.o : tagEventsNew.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

tagEventsNew : tagEventsNew.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

Fit_Daniel.o : Fit_Daniel.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

Fit_Daniel : Fit_Daniel.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

plotMuons.o : plotMuons.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

plotMuons : plotMuons.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skimZs.o : skimZs.C skimZs.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimZs : skimZs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitD0withZ.o : fitD0withZ.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitD0withZ : fitD0withZ.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYields.o : getTotalCharmYields.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYields : getTotalCharmYields.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYieldsMC.o : getTotalCharmYieldsMC.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYieldsMC : getTotalCharmYieldsMC.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYieldsDMC.o : getTotalCharmYieldsDMC.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYieldsDMC : getTotalCharmYieldsDMC.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYieldsJMC.o : getTotalCharmYieldsJMC.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYieldsJMC : getTotalCharmYieldsJMC.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYieldsWeighted.o : getTotalCharmYieldsWeighted.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYieldsWeighted : getTotalCharmYieldsWeighted.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYieldsJMCWeighted.o : getTotalCharmYieldsJMCWeighted.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYieldsJMCWeighted : getTotalCharmYieldsJMCWeighted.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeNewCalibTuples.o : makeNewCalibTuples.C makeNewCalibTuples.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeNewCalibTuples : makeNewCalibTuples.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeNewDTuples.o : makeNewDTuples.C makeNewDTuples.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeNewDTuples : makeNewDTuples.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getTotalCharmYieldsNew.o : getTotalCharmYieldsNew.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getTotalCharmYieldsNew : getTotalCharmYieldsNew.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeDDRPlots.o : makeDDRPlots.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeDDRPlots : makeDDRPlots.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitmcorn.o : fitmcorn.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitmcorn : fitmcorn.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

countGenXc.o : countGenXc.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

countGenXc : countGenXc.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitCharmNew.o : fitCharmNew.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitCharmNew : fitCharmNew.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

jpsiEff.o : jpsiEff.C jpsiEff.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

jpsiEff : jpsiEff.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getDPtDist.o : getDPtDist.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getDPtDist : getDPtDist.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

EffFromSV.o : EffFromSV.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

EffFromSV : EffFromSV.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

fitCharmFromSV.o : fitCharmFromSV.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitCharmFromSV : fitCharmFromSV.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeNewSVD0Effs.o : makeNewSVD0Effs.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeNewSVD0Effs : makeNewSVD0Effs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeNewSVD0EffsAdaptiveBinning.o : makeNewSVD0EffsAdaptiveBinning.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeNewSVD0EffsAdaptiveBinning : makeNewSVD0EffsAdaptiveBinning.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getMCD0StatsNew.o : getMCD0StatsNew.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getMCD0StatsNew : getMCD0StatsNew.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

plotD0Mass.o : plotD0Mass.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

plotD0Mass : plotD0Mass.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeNewD0Effs.o : makeNewD0Effs.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeNewD0Effs : makeNewD0Effs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeCombShape.o : makeCombShape.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeCombShape : makeCombShape.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
