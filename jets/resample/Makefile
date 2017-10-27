CXX = `root-config --cxx`

#`root-config --cxx --cflags` -o get_moments_B2Kstll get_moments_B2Kstll.C `root-config --glibs`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	resampleTrackIPs checkDaVinciSVs ipCorrelations

resampleTrackIPs.o : resampleTrackIPs.C resampleTrackIPs.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

resampleTrackIPs : resampleTrackIPs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

checkDaVinciSVs.o : checkDaVinciSVs.C checkDaVinciSVs.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

checkDaVinciSVs : checkDaVinciSVs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

ipCorrelations.o : ipCorrelations.C ipCorrelations.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

ipCorrelations : ipCorrelations.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
