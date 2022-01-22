CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer -lTMVA
INCPATH   = -I.

all : 	TreeSkimmer SimSkimmer

TreeSkimmer.o : TreeSkimmer.C TreeSkimmer.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

TreeSkimmer : TreeSkimmer.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

SimSkimmer.o : SimSkimmer.C SimSkimmer.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

SimSkimmer : SimSkimmer.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

EffSkimmer.o : EffSkimmer.C EffSkimmer.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

EffSkimmer : EffSkimmer.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

applyTMVA.o : applyTMVA.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

applyTMVA : applyTMVA.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
