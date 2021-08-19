#ifndef SVFITTER
#define SVFITTER

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TH1.h"
#include "TString.h"

class RooAbsData;
class RooAbsPdf;
class RooRealVar;

struct SVFitterOptions {
	int nMCorBins{94};
	double minMCor{600};
	double maxMCor{10000};
	int nNTrkBins{3};
	bool rerunTemplates{false};
	bool usePtBinnedTemplates{false};
	bool useBinnedTemplates{false};
	bool smearMCorShapes{false};
	bool additiveMCorCorrectionFactor{false};
	bool lightYieldFloat{false};
	double lightYieldScale{1.};
	bool correctBackTagEff{false};
	std::string useBackTagEffFromFile{""};
	bool runToyFits{false};
};

class SVFitter {
	public:
		SVFitter(TString name, TH1D* ptBins=0, TH1D* yBins=0) 
			: _name(name),
			_lightIsWeighted(false), _charmIsWeighted(true), _beautyIsWeighted(true),
			_inputsSet(false), _ptBins(ptBins), _yBins(yBins),
			_nmbins(18), _mmin(500.), _mmax(5000.), _ntbins(3) {}

		void setInputs(TString light, TString charm, TString beauty, TString data, bool lightIsMC=false, bool charmIsMC=true, bool beautyIsMC=true, bool dataIsMC=false, TString training="", TString back="");
		void setInputWeightings(bool light, bool charm, bool beauty);
		void setOptions(SVFitterOptions& options);
		void setSVBinning(int nmbins, double mmin, double mmax, int ntbins);
		void setRerunTemplates(bool rerun=true) {_recreateInputs = rerun;}
		void setUsePtBinnedTemplates(bool use=true) {_usePtBinnedTemplates = use;}
		void makeSVFitHists();
		void makeSVFitHists(int which);

		bool fit(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT=0u, uint binY=0u, TString enhancedTemplatesFile="");
		bool enhanceTemplatesAndFit(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT=0u, uint binY=0u, uint nPass=2);

		bool fitSVComb(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);
		bool fitSVSim(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u, uint enhance=0u, TString nameMod="");
		bool fitSV2D(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);
		bool fitSV2DUB(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);


		TString templateFileName(TString which);

	private:
		//calculate correction to named component(s) to match PDF to data
		//return the size of the largest deviation from unity (used as a criterion in iterative correction)
		double calcCorrection(RooRealVar* var, RooAbsData* ds, RooAbsPdf* pdf, TString toCorrect, TString plotName);
		void cachePDFs(RooRealVar* var, RooAbsData* ds, RooAbsPdf* pdf, TString toCorrect, TF1* corrFunc=0);

		void makeMassPlots(RooAbsData* ds, RooAbsPdf* pdf, TString name, int binPT, int binY, int minNBin, int maxNBin);

		TString _name;

		TString _lightInputFile;
		TString _charmInputFile;
		TString _beautyInputFile;
		TString _dataInputFile;
		TString _trainingInputFile;
		TString _backwardsInputFile;

		bool _lightIsWeighted;
		bool _charmIsWeighted;
		bool _beautyIsWeighted;

		bool _lightIsMC{false};
		bool _charmIsMC{true};
		bool _beautyIsMC{true};
		bool _dataIsMC{false};

		bool _inputsSet;

		TH1D* _ptBins;
		TH1D* _yBins;

		int _nmbins;
		double _mmin;
		double _mmax;
		int _ntbins;

		bool _recreateInputs{false};
		bool _usePtBinnedTemplates{false};
		bool _useBinnedTemplates{false};
		bool _smearMCorShapes{false};
		bool _additiveMCorCorrectionFactor{false};

		bool _lightYieldFloat{false};
		double _lightYieldScale{1.};
		bool _correctBackTagEff{false};
		TString _useBackTagEffFromFile{""};

		bool _runToyFits{false};

		bool _finalRun{false};
		int _nToys{100};

		//std::vector<std::vector<double>> _cShiftVals, _bShiftVals;
		std::vector<std::map<std::string,double>*> _fracYieldsCache;
		std::vector<std::map<std::string,double>*> _smearCache;
		std::vector<std::map<std::string,TH1*>*> _pdfCache;

		TString loadedTemplateFile;

		//bool _fixEnhancements{false};

};
#endif
