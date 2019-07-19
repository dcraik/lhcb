CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	skim skimNew lumi

skim.o : skim.C skim.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim : skim.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skimNew.o : skimNew.C skimNew.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimNew : skimNew.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

lumi.o : lumi.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

lumi : lumi.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
