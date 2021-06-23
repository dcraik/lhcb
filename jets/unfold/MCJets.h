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
			jetRecoDp4,
			jetRecoDp5,
			jetRecoDs4,
			jetRecoDs5,
			jetRecoLc4,
			jetRecoLc5,
			jetRecoSV4,
			jetRecoSV5,
			jetTrueD04,
			jetTrueD05,
			jetTrueDp4,
			jetTrueDp5
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
			  _dptmin(5000.), _dptmax(-1.),
			  _ybins(0), _recreateInputs(false) {}

		void setInputs(TString light, TString charm, TString beauty, TString svData="", TString d0Data="", TString dpData="");
		void setInputTruePtWeights(jetType type, TH1D* truePtWeights);
		void setExtraPtWeights(jetType type, TH1D* recoPtWeights);
		void useInputPtWeightsFromTree(bool useInputPtWeights=true) { _useInputPtWeights = useInputPtWeights;}
		void setDPtRange(double minpt, double maxpt);
		void setYBins(TH1D* ybins) { _ybins = ybins; }

		void setPtCorrFactorsFile(TString ptCorrFactorsFile) {_ptCorrFactorsFile = ptCorrFactorsFile;}
		void setRerunWeights(bool rerun=true) {_recreateInputs = rerun;}

		void setUnfoldingMethod(UnfoldMethod unfoldingMethod, double regPar);
		void setUnfoldingMethod(TString unfoldingMethod, double regPar);
		void setUnfoldingScaleSmear(double scaleFactor=1., double smearFactor=0.);

		double getPtCorrFactor(jetType type, double ptMin, double ptMax);

		bool weightMC(jetType type, TH2D* effTrueHist=0, TH2D* effRecoHist=0);

		TH1D* unfold(TH1D* input, jetType type, uint ybin=0);
		TH2D* unfold(TH2D* input, jetType type, bool useYBins=false);

		//TString setupTrainTestSamples(int sample);
		bool getTruth(TString file, std::map<std::string,TH1D*> hists, bool useTruePT=false);
		bool getTruth(TString file, TH1D* true4, TH1D* true5, TH1D* trueD04, TH1D* trueD05, TH1D* trueD0Sel4, TH1D* trueD0Sel5, TH1D* trueSV4, TH1D* trueSV5, bool useTruePT=false);
		bool getTruth(TString file, TH2D* true0, TH2D* true4, TH2D* true5, TH2D* trueSV4, TH2D* trueSV5, bool useTruePT=false);


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
		TString _dpDataInput;

		TString _resampledName;

		bool _inputsSet;
		bool _mcSampled;

		UnfoldMethod _unfoldingMethod{unfoldBayes};
		double _regPar{4};
		double _unfoldingScaleFactor{1.};
		double _unfoldingSmearFactor{0.};

		std::map<jetType,TH1D*> _truePtWeights;
		std::map<jetType,TH1D*> _extraPtWeights;
		std::map<jetType, std::map<uint,RooUnfoldResponse*> > _unfoldings;
		bool _useInputPtWeights;

		double _dptmin;
		double _dptmax;

		TH1D* _ybins;

		bool _recreateInputs;

		TString _ptCorrFactorsFile{""};
};

#endif
