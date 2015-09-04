CC  = c++
CXX = c++

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER
ROOTFLAGS = `root-config --cflags --libs` `gsl-config --cflags` -lHtml -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
INCPATH   = -Isrc/

all : 	objects fit fitEff

clean :
	rm -f bin/*.o

objects : 	bin/RooB2Kll.o bin/fit.o

bin/RooB2Kll.o : src/RooB2Kll.cxx src/RooB2Kll.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

bin/fit.o : src/fit.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fit : bin/fit.o bin/RooB2Kll.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^

bin/fitEff.o : src/fitEff.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

fitEff : bin/fitEff.o bin/RooB2Kll.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^
