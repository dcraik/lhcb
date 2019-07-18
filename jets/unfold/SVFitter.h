#ifndef SVFITTER
#define SVFITTER

#include <iostream>

#include "TH1.h"
#include "TString.h"

class SVFitter {
	public:
		SVFitter(TString name, TH1D* ptBins=0, TH1D* yBins=0) 
			: _name(name),
			_lightIsWeighted(false), _charmIsWeighted(true), _beautyIsWeighted(true),
			_inputsSet(false), _ptBins(ptBins), _yBins(yBins),
			_nmbins(18), _mmin(500.), _mmax(5000.), _ntbins(3) {}

		void setInputs(TString light, TString charm, TString beauty, TString data);
		void setInputWeightings(bool light, bool charm, bool beauty);
		void setSVBinning(int nmbins, double mmin, double mmax, int ntbins);
		void makeSVFitHists(int which);
		bool fitSV(double& NB, double& eB, double& NC, double& eC, double& NQ, double& eQ, uint binPT, uint binY=1u);

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

		bool _inputsSet;

		TH1D* _ptBins;
		TH1D* _yBins;

		int _nmbins;
		double _mmin;
		double _mmax;
		int _ntbins;
};
#endif
