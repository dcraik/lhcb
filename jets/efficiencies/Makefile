CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	makeNewD0Effs addEfficiencies

makeNewD0Effs.o : makeNewD0Effs.C makeNewD0Effs.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeNewD0Effs : makeNewD0Effs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

addEfficiencies.o : addEfficiencies.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

addEfficiencies : addEfficiencies.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getMCD0StatsNew.o : getMCD0StatsNew.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getMCD0StatsNew : getMCD0StatsNew.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
