#ifndef MCJETS_H
#define MCJETS_H

#include <map>

#include "TString.h"

class TH1D;
class TH2D;

class RooUnfoldResponse;

class MCJets {
	public:
		enum jetType {
			jetAll0,
			jetAll4,
			jetAll5,
			jetRecoD04,
			jetRecoD05,
			jetRecoSV4,
			jetRecoSV5
		};

		enum UnfoldMethod {
			unfoldInvert,
			unfoldBayes,
			unfoldIDS,
			unfoldSVD,
			unfoldTUnfold
		};

		MCJets(TString name)
			: _name(name), _inputsSet(false), _mcSampled(false),
			  _d0ptmin(5000.), _d0ptmax(-1.),
			  _ybins(0), _recreateInputs(false), _useFastCorrFactors(false) {}

		void setInputs(TString light, TString charm, TString beauty, TString svData="", TString d0Data="");
		void setInputTruePtWeights(jetType type, TH1D* truePtWeights);
		void setD0PtRange(double minpt, double maxpt);
		void setYBins(TH1D* ybins) { _ybins = ybins; }

		void setFastCorrFactors(bool useFastCorrFactors) {_useFastCorrFactors = useFastCorrFactors;}
		void setRerunWeights(bool rerun=true) {_recreateInputs = rerun;}

		void setUnfoldingMethod(UnfoldMethod unfoldingMethod, double regPar);
		void setUnfoldingMethod(TString unfoldingMethod, double regPar);

		double getPtCorrFactor(jetType type, double ptMin, double ptMax);

		bool weightMC(jetType type, TH2D* effTrueHist=0, TH2D* effRecoHist=0);

		TH1D* unfold(TH1D* input, jetType type, uint ybin=0);
		TH2D* unfold(TH2D* input, jetType type, bool useYBins=false);

		//TString setupTrainTestSamples(int sample);
		bool getTruth(TString file, TH1D* true4, TH1D* true5, TH1D* trueD04, TH1D* trueD05, TH1D* trueD0Sel4, TH1D* trueD0Sel5, TH1D* trueSV4, TH1D* trueSV5, bool useTruePT=false);

		const TString outputName(jetType type);
		const TString typeName(jetType type);
		const TString resampledName() { return _resampledName; }
	private:
		RooUnfoldResponse* trainUnfold(TH1D* binning, jetType type, uint ybin=0);

		TString _name;

		TString _lightInput;
		TString _charmInput;
		TString _beautyInput;
		TString _svDataInput;
		TString _d0DataInput;

		TString _resampledName;

		bool _inputsSet;
		bool _mcSampled;

		UnfoldMethod _unfoldingMethod{unfoldBayes};
		double _regPar{4};

		std::map<jetType,TH1D*> _truePtWeights;
		std::map<jetType, std::map<uint,RooUnfoldResponse*> > _unfoldings;

		double _d0ptmin;
		double _d0ptmax;

		TH1D* _ybins;

		bool _recreateInputs;

		bool _useFastCorrFactors;
};

#endif
