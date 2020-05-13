#ifndef SVFITTER
#define SVFITTER

#include <iostream>

#include "TH1.h"
#include "TString.h"

struct SVFitterOptions {
	int nMCorBins{94};
	double minMCor{600};
	double maxMCor{10000};
	int nNTrkBins{3};
	bool rerunTemplates{false};
	bool usePtBinnedTemplates{false};
	bool useBinnedTemplates{false};
};

class SVFitter {
	public:
		SVFitter(TString name, TH1D* ptBins=0, TH1D* yBins=0) 
			: _name(name),
			_lightIsWeighted(false), _charmIsWeighted(true), _beautyIsWeighted(true),
			_inputsSet(false), _ptBins(ptBins), _yBins(yBins),
			_nmbins(18), _mmin(500.), _mmax(5000.), _ntbins(3) {}

		void setInputs(TString light, TString charm, TString beauty, TString data, bool lightIsMC=false, bool charmIsMC=true, bool beautyIsMC=true);
		void setInputWeightings(bool light, bool charm, bool beauty);
		void setOptions(SVFitterOptions& options);
		void setSVBinning(int nmbins, double mmin, double mmax, int ntbins);
		void setRerunTemplates(bool rerun=true) {_recreateInputs = rerun;}
		void setUsePtBinnedTemplates(bool use=true) {_usePtBinnedTemplates = use;}
		void makeSVFitHists(int which);
		bool fitSV(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);
		bool fitSVSim(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);
		bool fitSV2D(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);
		bool fitSV2DUB(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);

		TString templateFileName(TString which);

	private:
		TString _name;

		TString _lightInputFile;
		TString _charmInputFile;
		TString _beautyInputFile;
		TString _dataInputFile;

		bool _lightIsWeighted;
		bool _charmIsWeighted;
		bool _beautyIsWeighted;

		bool _lightIsMC{false};
		bool _charmIsMC{true};
		bool _beautyIsMC{true};

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
};
#endif
