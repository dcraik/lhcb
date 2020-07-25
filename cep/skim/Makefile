CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	skim skimMC getLumi getLumiFromZ addPtEffWeight

skim.o : skim.C skim.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim : skim.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skimMC.o : skimMC.C skimMC.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimMC : skimMC.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getLumi.o : getLumi.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getLumi : getLumi.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

getLumiFromZ.o : getLumiFromZ.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

getLumiFromZ : getLumiFromZ.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

addPtEffWeight.o : addPtEffWeight.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

addPtEffWeight : addPtEffWeight.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
