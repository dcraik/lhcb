CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I.

all : 	skim_for_yandex	skim_for_yandex_data fit fitD0 fitD fitDs fitLc fit2d

skim_for_yandex.o : skim_for_yandex.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex : skim_for_yandex.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

skim_for_yandex_data.o : skim_for_yandex_data.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

skim_for_yandex_data : skim_for_yandex_data.o
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

