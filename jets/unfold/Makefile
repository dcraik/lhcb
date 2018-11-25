CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I. -I./RooUnfold/src

all : 	unfoldJetPt RooUnfoldD0Mass

RooPromptShape.o : RooPromptShape.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

unfoldJetPt.o : unfoldJetPt.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

unfoldJetPt : unfoldJetPt.o RooUnfold/libRooUnfold.so
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

testUnfolding.o : testUnfolding.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testUnfolding : testUnfolding.o RooUnfold/libRooUnfold.so
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

RooUnfoldD0Mass.o : RooUnfoldD0Mass.cxx 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

RooUnfoldD0Mass : RooUnfoldD0Mass.o RooUnfold/libRooUnfold.so
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

correctAndUnfold.o : correctAndUnfold.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

correctAndUnfold : correctAndUnfold.o RooUnfold/libRooUnfold.so
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

fitCorrectAndUnfold.o : fitCorrectAndUnfold.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

fitCorrectAndUnfold : fitCorrectAndUnfold.o RooUnfold/libRooUnfold.so RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

fitSVAndUnfold.o : fitSVAndUnfold.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

fitSVAndUnfold : fitSVAndUnfold.o RooUnfold/libRooUnfold.so RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

getTaggingEffs.o : getTaggingEffs.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

getTaggingEffs : getTaggingEffs.o RooUnfold/libRooUnfold.so RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

massFitExample.o : massFitExample.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

massFitExample : massFitExample.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so
