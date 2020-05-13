#pragma once

#include <string>
#include <vector>

#include "TString.h"

class TH1D;
class TTree;

class RooAbsData;
class RooAbsPdf;
class RooArgList;
class RooRealVar;
class RooWorkspace;

struct SimDFitter {
	enum class fitType {
		fitMass,
		fitIP,
		fit2D
	};
	enum class fitSet {
		fitJustCharm,
		fitJustBeauty,
		fitJustSignal,
		fitJustSideband,
		fitData,
		fitToys
	};
	enum class truthMatchType {
		truthMatchOff,
		truthMatchPrompt=4,
		truthMatchDispl=5
	};

	SimDFitter(TString name, TH1D* ptBins, TH1D* yBins=0) 
		: _name(name) { initJetBinning(ptBins, yBins);}

	void setInputs(TString data, TString charm, TString beauty, TString eff, TString acc, bool isMC=false);
	void setDPtRange(double ptmin, double ptmax);
	TString const dFileName();

	inline void setSampleType(fitSet sample) { _sample = sample; }
	inline void setTruthMatchType(truthMatchType type) { _truthMatchData = type; }

	bool addEffs();
	bool addEffs(int flavour);
	bool testEffs(int flavour);
	bool fitD(double& yield4, double& error4, double& yield5, double& error5, uint binPT, uint binY=0u, uint effType=4u);
	double getAveWeight(double flavour);

	void setRerunEffs(bool rerun=true) {_recreateInputs = rerun;}
	void skipSumW2Fits(bool skip=true) {_skipSumW2Fits = skip;}

private:
	void initJetBinning(TH1D* ptBins, TH1D* yBins);

	bool init();
	//bool makeInput(int flav);
	bool loadDatasets();
	bool setupDPtBins(TTree* t);
	bool setupModel();
	void initPars();

	bool fit(fitType whichFit, fitSet whichSet, uint ptBinMin, uint ptBinMax, unsigned int nToys=0);
	void makeMassPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin);
	void makeIPPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin);
	void make2DPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin);

	bool doFits();

	std::vector<double> _ptBinsD;//{0.,7e3,10e3,15e3,20e3,100e3};
	bool _usePtFracBins{true};
	uint _nDPtBins{5u};

	std::vector<double> _ptBins{0.,100e3};
	std::vector<double> _yBins{0.,10.};

	double valMassMean{1864.};
	double valMassWidth{10.};
	double valMassRatio{1.5};
	double valMassAlpha{2.};
	double valMassN{2.};
	double valMassFrac{0.7};
	double valBkgFrac{0.1};

	double valPromptMean{0.7};
	double valPromptWidth{1.1};
	double valPromptAsym{-0.25};
	double valPromptRhoL{1.3};
	double valPromptRhoR{1.7};
	double valDisplMean{4.};
	double valDisplWidth{3.};
	double valDisplAsym{-0.29};
	double valDisplRhoL{5.};
	double valDisplRhoR{5.};

	double _yield4{0.};
	double _error4{0.};
	double _yield5{0.};
	double _error5{0.};

	RooWorkspace* ws{0};

	double nEntries{0};

	TString _dataFile;
	TString _combFile;
	TString _charmFile;
	TString _beautyFile;
	TString _effFile;
	TString _accFile;

	fitSet _sample{fitSet::fitData};
	bool _dataIsMC{false};
	truthMatchType _truthMatchData{truthMatchType::truthMatchOff};

	bool _recreateInputs{false};

	bool _inputsSet{false};
	bool _effsSet{false};
	bool _initialised{false};
	bool _mcFitsDone{false};

	TString _name;

	double _dptmin{5000.};
	double _dptmax{-1.};

	//parameters set by fitD to be used in private doFits function
	int _flavour;
	double _jptmin;
	double _jptmax;
	double _ymin;
	double _ymax;
	bool _useEffs;
	bool _skipSumW2Fits{false};
};
