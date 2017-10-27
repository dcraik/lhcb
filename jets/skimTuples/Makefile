CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	skim_for_yandex	skim_for_yandex_data skim_for_yandex_data_z0 skim_for_yandex_data_light skim_for_yandex_DMC skimD0s skimD0sFromB skimCharm charmEfficiencies weightCharmEffs fit fitD0 fitD fitDs fitLc fit2d tagEvents skimZs

skim_for_yandex.o : skim_for_yandex.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex : skim_for_yandex.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data.o : skim_for_yandex_data.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data : skim_for_yandex_data.o
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

RooPromptShape.o : RooPromptShape.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

tagEvents.o : tagEvents.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

tagEvents : tagEvents.o
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
