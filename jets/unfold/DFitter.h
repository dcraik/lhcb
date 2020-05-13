#ifndef DFITTER
#define DFITTER

#include <iostream>

#include "TH1.h"
#include "TString.h"

class DFitter {
	public:
		enum fitType{
			fit1D1D,
			fit2D,
			fitSidebandSub,
			fitSBS1D
		};
		enum sampleType{
			mixedSample,
			charmSample,
			beautySample
		};

		DFitter(TString name, TH1D* ptBins, TH1D* yBins=0) 
			: _name(name), _type(fitSBS1D), _sample(mixedSample),
			_useSimpleEffs(false), _dataIsMC(false),
			_useRhoZCorr(true),
			_inputsSet(false), _effsSet(false),
			_ptBins(ptBins), _yBins(yBins),
			_dptmin(5000.), _dptmax(-1.),
			_flavour(4), _jptmin(0.), _jptmax(1e5),
			_ymin(0.), _ymax(10.) {}

		void setInputs(TString data, TString eff, TString acc, bool simple=false, bool isMC=false, bool useRhoZCorr=true);
		void setDPtRange(double ptmin, double ptmax);
		inline void setFitType(fitType type) { _type = type; }
		inline void setSampleType(sampleType sample) { _sample = sample; }

		bool addEffs();
		bool testEffs(int flavour);
		bool fitD(int flavour, double& yield, double& error, uint binPT, uint binY=1u, bool useEffs=true);

		TString const dFileName();

		TString const typeName(fitType type);

	private:
		bool addEffsFull();
		bool addEffsSimple();

		bool testEffsFull();
		bool testEffsSimple();

		bool fitD_1D1D(double& yield, double& error);
		bool fitD_SBS(double& yield, double& error);
		bool fitD_SBS1D(double& yield, double& error);
		bool fitD_2D(double& yield, double& error);

		TString _name;

		fitType _type;
		sampleType _sample;

		TString _dataFile;
		TString _effFile;
		TString _accFile;

		bool _useSimpleEffs;
		bool _dataIsMC;
		bool _useRhoZCorr;

		bool _inputsSet;
		bool _effsSet;

		TH1D* _ptBins;
		TH1D* _yBins;

		double _dptmin;
		double _dptmax;

		//parameters set by fitD to be used in private fitD_ functions
		int _flavour;
		double _jptmin;
		double _jptmax;
		double _ymin;
		double _ymax;
		bool _useEffs;
};
#endif
