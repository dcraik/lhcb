CC  = c++
CXX = c++

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER
ROOTFLAGS = `root-config --cflags --libs` `gsl-config --cflags` -lHtml -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
INCPATH   = -Isrc/

all : 	objects gen

clean :
	rm -f bin/*.o

objects : 	bin/RooB2Kll.o bin/gen.o

bin/RooB2Kll.o : src/RooB2Kll.cxx src/RooB2Kll.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

bin/gen.o : src/gen.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

gen : bin/gen.o bin/RooB2Kll.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^
