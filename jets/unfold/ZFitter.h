#ifndef ZFITTER
#define ZFITTER

#include <iostream>

#include "TH1.h"
#include "TString.h"

class ZFitter {
	public:
		ZFitter(TString name, TH1D* ptBins, TH1D* yBins) 
			: _name(name),
			_inputsSet(false), _fixed(false),
			_ptBins(ptBins), _yBins(yBins),
			_jptmin(0.), _jptmax(1e5),
			_ymin(0.), _ymax(10.),
			_ratio(1.), _frac(0.5), _alpha1(2.), _alpha2(-2.) {}

		void setInputs(TString data);

		bool fit(double& yield, double& error, uint binPT=0, uint binY=0);
		void fixShape(bool fixed=true) { _fixed=fixed; }

	private:
		TString _name;

		TString _dataFile;

		bool _inputsSet;
		bool _fixed;

		TH1D* _ptBins;
		TH1D* _yBins;

		//parameters set by fit
		double _jptmin;
		double _jptmax;
		double _ymin;
		double _ymax;

		//shape parameters
		double _ratio;
		double _frac;
		double _alpha1;
		double _alpha2;
};
#endif
