CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I. -I./RooUnfold/src

all : 	unfoldJetPt RooUnfoldD0Mass

cloneHists.o : cloneHists.C cloneHists.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

outputFunctions.o : outputFunctions.C outputFunctions.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

SVFitter.o : SVFitter.C SVFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

DFitter.o : DFitter.C DFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

ZFitter.o : ZFitter.C ZFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

MCJets.o : MCJets.C MCJets.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

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

getTaggingEffs : getTaggingEffs.o RooUnfold/libRooUnfold.so RooPromptShape.o outputFunctions.o cloneHists.o SVFitter.o DFitter.o MCJets.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

massFitExample.o : massFitExample.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

massFitExample : massFitExample.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

fitZJetSVs.o : fitZJetSVs.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

fitZJetSVs : fitZJetSVs.o outputFunctions.o cloneHists.o SVFitter.o ZFitter.o MCJets.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

getCharmBkgFracs.o : getCharmBkgFracs.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

getCharmBkgFracs : getCharmBkgFracs.o outputFunctions.o cloneHists.o SVFitter.o MCJets.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

fitSVToys.o : fitSVToys.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

fitSVToys : fitSVToys.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

quickSVToys.o : quickSVToys.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

quickSVToys : quickSVToys.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

testSimPtWeighting.o : testSimPtWeighting.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testSimPtWeighting : testSimPtWeighting.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
