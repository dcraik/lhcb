#ifndef DOFITS_H
#define DOFITS_H

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooLinkedList.h"

#include "TString.h"

#include <vector>

class TChain;
class TFile;
class TTree;

class Fitter {
	public:
		enum fitType {
			fitSim,
			fitpp,
			fitppNoPrevEtCut,
			fitpp5,
			fitpPb,
			fitPbp,
			fitPbPb
		};

		static const int optNone{0x00};
		static const int optHiPTSq{0x01};
		static const int optHiHRC{0x02};
		static const int optHiNTrk{0x04};

		Fitter() = default;

		void setup();
		void run();

		TString const typeName(fitType type);
		TString const optsName(uint options);
		TString const optsCut(uint options);

	private:
		void fitEff();
		void fitMass(fitType which, uint options=optNone, bool comb=false);//bool hiPTSq=false, bool hiHRC=false, bool hiNTrk=false, bool comb=false);
		void fitPtSq(fitType which, uint options=optNone);//bool hiPTSq=false, bool hiHRC=false, bool hiNTrk=false);

		//helper function to add constraints to a fit
		bool addGaussConstraint(RooArgSet& pdfs, RooLinkedList& fitCmds, RooAbsReal& var, double mean, double sigma);

		TTree* getReducedTree(fitType which, uint options, bool massFit, bool comb);//, bool hiPTSq, bool hiHRC, bool hiNTrk, bool comb);

		//mass fit parameters
		double _valMu{1020.}, _valSigma{1.5};
		//double _valA1{2.}, _valA2{-2.}, _valN1{2.}, _valN2{2.}, _valR{1.1}, _valF{0.5};
		double _valC1{-0.03}, _valC2{980.};

		//efficiency parameters
		//double _x0Val{0.331}, _x0Error{0.018};
		//double _x1Val{3.50},  _x1Error{0.14};
		//double _gLVal{3.02},  _gLError{0.19};
		//double _gRVal{0.68},  _gRError{0.04};
		double _x0Val{0.314}, _x0Error{0.010};
		double _x1Val{3.31},  _x1Error{0.07};
		double _gLVal{3.08},  _gLError{0.11};
		double _gRVal{0.709}, _gRError{0.026};

		std::vector<double> _splineXVals;
		std::vector<double> _splineYVals;

		//pT^2 fit parameters
		double _bkgYield{0.}, _bkgError{0.}, _aDis{-1.}, _aDisError{0.};

		bool _doBinnedFit{true};
		bool _useSplineEff{true};
		bool _floatEffFunc{false};

		//datasets
		TChain *_ppData{0}, *_pp5Data{0}, *_ppbData{0}, *_pbpData{0}, *_pbpbData{0}, *_simData{0}, *_ppCmbData{0}, *_pp5CmbData{0}, *_ppbCmbData{0}, *_pbpCmbData{0}, *_pbpbCmbData{0};

		//file for reduced trees used in fits
		TFile *_fitTreeFile{0};
};

#endif
