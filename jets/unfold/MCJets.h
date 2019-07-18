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

		MCJets(TString name)
			: _name(name), _inputsSet(false),
			  _d0ptmin(5000.), _d0ptmax(-1.),
			  _ybins(0) {}

		void setInputs(TString light, TString charm, TString beauty, TString svData="", TString d0Data="");
		void setInputTruePtWeights(jetType type, TH1D* truePtWeights);
		void setD0PtRange(double minpt, double maxpt);
		void setYBins(TH1D* ybins) { _ybins = ybins; }

		bool weightMC(jetType type, TH2D* effHist=0);

		TH1D* unfold(TH1D* input, jetType type, uint ybin=0);
		TH2D* unfold(TH2D* input, jetType type, bool useYBins=false);

		const TString outputName(jetType type);
		const TString typeName(jetType type);
	private:
		RooUnfoldResponse* trainUnfold(TH1D* binning, jetType type, uint ybin=0);

		TString _name;

		TString _lightInput;
		TString _charmInput;
		TString _beautyInput;
		TString _svDataInput;
		TString _d0DataInput;

		bool _inputsSet;

		std::map<jetType,TH1D*> _truePtWeights;
		std::map<jetType, std::map<uint,RooUnfoldResponse*> > _unfoldings;

		double _d0ptmin;
		double _d0ptmax;

		TH1D* _ybins;
};

#endif
