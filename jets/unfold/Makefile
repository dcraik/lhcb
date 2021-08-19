CXX = `root-config --cxx`

CXXFLAGS  = -O2 -Wall -Wno-deprecated -D__ROOFIT_NOBANNER -std=c++11 -g
ROOTFLAGS = `root-config --cflags` # --libs` `gsl-config --cflags` -lThread -lMinuit -lMathCore -lMinuit2 -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer 
ROOTLIBS  = `root-config --glibs` -lRooFit -lThread -lMinuit -lMathCore -lRooStats -lRooFit -lRooFitCore -lFoam -lTreePlayer
INCPATH   = -I. -I./RooUnfold/src

all : 	fitZJetSVs getTaggingEffs

clean :
	rm -rf *.o fitZJetSVs getTaggingEffs

cloneHists.o : cloneHists.C cloneHists.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

AdaptBin.o : AdaptBin.cc AdaptBin.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

outputFunctions.o : outputFunctions.C outputFunctions.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

SVFitter.o : SVFitter.C SVFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

DFitter.o : DFitter.C DFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

SimDFitter.o : SimDFitter.C SimDFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

ZFitter.o : ZFitter.C ZFitter.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

MCJets.o : MCJets.C MCJets.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

DatasetManager.o : DatasetManager.C DatasetManager.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testDM : DatasetManager.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

RooPromptShape.o : RooPromptShape.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

RooCatKeysPdf.o : RooCatKeysPdf.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

RooMultiKeysPdf.o : RooMultiKeysPdf.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

RooModulatedHistPdf.o : RooModulatedHistPdf.cxx
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

svTemplatePlots.o : svTemplatePlots.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

svTemplatePlots : svTemplatePlots.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

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

getTaggingEffs : getTaggingEffs.o RooUnfold/libRooUnfold.so RooModulatedHistPdf.o RooMultiKeysPdf.o RooPromptShape.o outputFunctions.o cloneHists.o SVFitter.o SimDFitter.o MCJets.o DatasetManager.o AdaptBin.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

testUnfoldingMC.o : testUnfoldingMC.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testUnfoldingMC : testUnfoldingMC.o RooUnfold/libRooUnfold.so outputFunctions.o cloneHists.o MCJets.o DatasetManager.o AdaptBin.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

testUnfoldingMC2.o : testUnfoldingMC2.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testUnfoldingMC2 : testUnfoldingMC2.o RooUnfold/libRooUnfold.so outputFunctions.o cloneHists.o MCJets.o DatasetManager.o AdaptBin.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

testUnfoldingMC3.o : testUnfoldingMC3.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testUnfoldingMC3 : testUnfoldingMC3.o RooUnfold/libRooUnfold.so outputFunctions.o cloneHists.o MCJets.o DatasetManager.o AdaptBin.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

testUnfoldingMC_compare.o : testUnfoldingMC_compare.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testUnfoldingMC_compare : testUnfoldingMC_compare.o RooUnfold/libRooUnfold.so outputFunctions.o cloneHists.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

testUnfoldingMC4.o : testUnfoldingMC4.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testUnfoldingMC4 : testUnfoldingMC4.o RooUnfold/libRooUnfold.so outputFunctions.o cloneHists.o MCJets.o DatasetManager.o AdaptBin.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

doDFits_DPtBins.o : doDFits_DPtBins.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

doDFits_DPtBins : doDFits_DPtBins.o RooUnfold/libRooUnfold.so RooModulatedHistPdf.o RooMultiKeysPdf.o RooPromptShape.o outputFunctions.o cloneHists.o SVFitter.o SimDFitter.o MCJets.o DatasetManager.o AdaptBin.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

d0FitTests.o : d0FitTests.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

d0FitTests : d0FitTests.o RooUnfold/libRooUnfold.so RooPromptShape.o outputFunctions.o cloneHists.o DFitter.o MCJets.o DatasetManager.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

massFitExample.o : massFitExample.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

massFitExample : massFitExample.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) RooUnfold/libRooUnfold.so

fitZJetSVs.o : fitZJetSVs.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

fitZJetSVs : fitZJetSVs.o outputFunctions.o cloneHists.o SVFitter.o ZFitter.o MCJets.o DatasetManager.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options RooUnfold/libRooUnfold.so

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

simJetTrigEffs.o : simJetTrigEffs.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

simJetTrigEffs : simJetTrigEffs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

simpleD0ToyStudy.o : simpleD0ToyStudy.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

simpleD0ToyStudy : simpleD0ToyStudy.o outputFunctions.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

twoStageD0ToyStudy.o : twoStageD0ToyStudy.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

twoStageD0ToyStudy : twoStageD0ToyStudy.o outputFunctions.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

simFitD0ToyStudy.o : simFitD0ToyStudy.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

simFitD0ToyStudy : simFitD0ToyStudy.o outputFunctions.o RooPromptShape.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

testRooCatKeysPdf.o : testRooCatKeysPdf.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

testRooCatKeysPdf : testRooCatKeysPdf.o outputFunctions.o RooCatKeysPdf.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeWeightSimTuple.o : makeWeightSimTuple.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

makeWeightSimTuple : makeWeightSimTuple.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeWeightSim5Tuple.o : makeWeightSim5Tuple.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

makeWeightSim5Tuple : makeWeightSim5Tuple.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeWeightSimTuple_step2.o : makeWeightSimTuple_step2.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

makeWeightSimTuple_step2 : makeWeightSimTuple_step2.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

makeWeightSim5Tuple_step2.o : makeWeightSim5Tuple_step2.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

makeWeightSim5Tuple_step2 : makeWeightSim5Tuple_step2.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

addWeightsToSim.o : addWeightsToSim.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

addWeightsToSim : addWeightsToSim.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

addWeightsToSim5.o : addWeightsToSim5.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

addWeightsToSim5 : addWeightsToSim5.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

SVFitSysts.o : SVFitSysts.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

SVFitSysts : SVFitSysts.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

calcJetRecoEffs.o : calcJetRecoEffs.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

calcJetRecoEffs : calcJetRecoEffs.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

effClosureTest.o : effClosureTest.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

effClosureTest : effClosureTest.o cloneHists.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)

effClosureTest_Dp.o : effClosureTest_Dp.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCPATH) -c $< -o $@

effClosureTest_Dp : effClosureTest_Dp.o cloneHists.o
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $^ $(ROOTLIBS)
