CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	skimTuples

velo.o : velo.cc velo.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimTuples.o : skimTuples.C skimTuples.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skimTuples : skimTuples.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeCombShape.o : makeCombShape.C 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

makeCombShape : makeCombShape.o velo.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
