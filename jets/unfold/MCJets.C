#include "MCJets.h"

#include <iostream>
#include <vector>

#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldDagostini.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"
#include "RooUnfoldInvert.h"

#include "cloneHists.h"
#include "outputFunctions.h"
#include "DatasetManager.h"

void MCJets::setInputs(TString light, TString charm, TString beauty, TString svData, TString d0Data, TString dpData) {
	_lightInput  = light;
	_charmInput  = charm;
	_beautyInput = beauty;
	_svDataInput = svData;
	_d0DataInput = d0Data;
	_dpDataInput = dpData;

	_inputsSet = true;
}

void MCJets::setInputTruePtWeights(jetType type, TH1D* truePtWeights) {
	_truePtWeights[type] = truePtWeights;
}

void MCJets::setExtraPtWeights(jetType type, TH1D* recoPtWeights) {
	_extraPtWeights[type] = recoPtWeights;
}

void MCJets::setDPtRange(double minpt, double maxpt) {
	_dptmin = minpt;
	_dptmax = maxpt;
}

void MCJets::setUnfoldingMethod(UnfoldMethod unfoldingMethod, double regPar=-1) {
	_unfoldingMethod = unfoldingMethod;
	_regPar = regPar;
}

void MCJets::setUnfoldingMethod(TString unfoldingMethod, double regPar=-1) {
	if(unfoldingMethod=="invert") {
		setUnfoldingMethod(unfoldInvert,regPar);
	} else if(unfoldingMethod=="bayes") {
		setUnfoldingMethod(unfoldBayes,regPar);
	} else if(unfoldingMethod=="ids") {
		setUnfoldingMethod(unfoldIDS,regPar);
	} else if(unfoldingMethod=="svd") {
		setUnfoldingMethod(unfoldSVD,regPar);
	} else if(unfoldingMethod=="tunfold") {
		setUnfoldingMethod(unfoldTUnfold,regPar);
	} else {
		std::cout << "Unknown unfolding method " << unfoldingMethod << std::endl;
	}
}

void MCJets::setUnfoldingScaleSmear(double scaleFactor, double smearFactor) {
		_unfoldingScaleFactor = scaleFactor;
		_unfoldingSmearFactor = smearFactor;
}

double MCJets::getPtCorrFactor(jetType type, double ptMin, double ptMax) {
	TString cut1;
	TString cut2;

	cut1 = "JetPT>"; cut1+=ptMin; cut1+=" && JetPT<"; cut1+=ptMax;
	cut2 = "JetPT>"; cut2+=ptMin; cut2+=" && JetPT<"; cut2+=ptMax;
	DatasetManager* dm = DatasetManager::getInstance();

	switch(type) {
		case jetAll0:
			std::cout << "getting pT correction factor for all jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetAll4:
			std::cout << "getting pT correction factor for all c-jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_charmInput);
			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000.";
			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0]";
			break;
		case jetAll5:
			std::cout << "getting pT correction factor for all b-jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000.";
			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0]";
			break;
		case jetRecoD04:
			std::cout << "getting pT correction factor for c->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetRecoDp4:
			std::cout << "getting pT correction factor for c->D+ in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetRecoDs4:
			std::cout << "getting pT correction factor for c->Ds in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetRecoLc4:
			std::cout << "getting pT correction factor for c->Lc in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetRecoSV4:
			std::cout << "getting pT correction factor for c->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_charmInput);
			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000. && SVM[0]";
			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && SVM[0]";
			break;
		case jetTrueD04:
			std::cout << "getting pT correction factor for c->all D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetTrueDp4:
			std::cout << "getting pT correction factor for c->all D+ in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetRecoD05:
			std::cout << "getting pT correction factor for b->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			//cut1+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEBPT>5000.";
			//cut2+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEDPT>5000.";
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetRecoDp5:
			std::cout << "getting pT correction factor for b->D+ in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==411 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==411 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetRecoDs5:
			std::cout << "getting pT correction factor for b->Ds in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==431 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==431 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetRecoLc5:
			std::cout << "getting pT correction factor for b->Lc in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==4122 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==4122 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetRecoSV5:
			std::cout << "getting pT correction factor for b->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000. && SVM[0]";
			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && SVM[0]";
			break;
		case jetTrueD05:
			std::cout << "getting pT correction factor for b->all D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetTrueDp5:
			std::cout << "getting pT correction factor for b->all D+ in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			dm->setDataset(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==411 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==411 && TRUEDTRK2IDX!=-1 && TRUEDTRK3IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		default:
			return 1.;
	}

	dm->reset();

	//Take mean from previous toy results to save time
	if(_ptCorrFactorsFile!="") {
		TFile* f = TFile::Open(_ptCorrFactorsFile);
		if(f) {
			TString tname = TString::Format("%s_%.0f_%.0f", typeName(type).Data(), ptMin, ptMax);
			TTree* t = static_cast<TTree*>(f->Get(tname));
			if(t) {
				TH1D htemp("htemp","",100,t->GetMinimum("corr"),t->GetMaximum("corr"));
				t->Draw("corr>>htemp");
				std::cout << "Prec-calculated orrection factor: " << htemp.GetMean() << std::endl;
				return htemp.GetMean();
			}
		}
	}

	//See if we've cached this result already
	TFile* f = TFile::Open(gSaveDir+"/ptCorrFactors.root","UPDATE");
	TString tname = TString::Format("%s_%.0f_%.0f", typeName(type).Data(), ptMin, ptMax);
	TTree* t = static_cast<TTree*>(f->Get(tname));
	double corr(1.), num(0.), denom(0.);
	if(!t) {
		//If the tree doesn't exist make it and calculate the correction factor
		t = new TTree(tname,"");

		t->Branch("corr",  &corr);
		t->Branch("num",   &num);
		t->Branch("denom", &denom);

		std::vector<double> counts = dm->getEntries(std::vector<TString>{cut1,cut2});
		corr  = counts[0]/counts[1];
		num   = counts[0];
		denom = counts[1];
		t->Fill();
		t->AutoSave();
	} else {
		t->SetBranchAddress("corr",  &corr);
		t->GetEntry(0);
	}

	std::cout << "Correction factor: " << corr << std::endl;
	f->Close();
	return corr;
}

//double MCJets::getPtCorrFactor(jetType type, double ptMin, double ptMax) {
//	TFile* f(0);
//	TTree* t(0);
//	TString cut1;
//	TString cut2;
//	double corr(1.);
//
//	cut1 = "JetPT>"; cut1+=ptMin; cut1+=" && JetPT<"; cut1+=ptMax;
//	cut2 = "JetPT>"; cut2+=ptMin; cut2+=" && JetPT<"; cut2+=ptMax;
//
//	switch(type) {
//		case jetAll0:
//			std::cout << "getting pT correction factor for all jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			return 1.;//no correction needed
//		case jetAll4:
//			std::cout << "getting pT correction factor for all c-jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			f = TFile::Open(_charmInput);
//			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000.";
//			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0]";
//			break;
//		case jetAll5:
//			std::cout << "getting pT correction factor for all b-jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			f = TFile::Open(_beautyInput);
//			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000.";
//			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0]";
//			break;
//		case jetRecoD04:
//			std::cout << "getting pT correction factor for c->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			return 1.;//no correction needed
//		case jetRecoSV4:
//			std::cout << "getting pT correction factor for c->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			f = TFile::Open(_charmInput);
//			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000. && SVM[0]";
//			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && SVM[0]";
//			break;
//		case jetRecoD05:
//			std::cout << "getting pT correction factor for b->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			f = TFile::Open(_beautyInput);
//			//cut1+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEBPT>5000.";
//			//cut2+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEDPT>5000.";
//			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
//			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
//			break;
//		case jetRecoSV5:
//			std::cout << "getting pT correction factor for b->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
//			f = TFile::Open(_beautyInput);
//			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000. && SVM[0]";
//			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && SVM[0]";
//			break;
//		default:
//			return 1.;
//	}
//
//	if(!f) return 1.;
//	t = dynamic_cast<TTree*>(f->Get("T"));
//	if(!t) return 1.;
//
//	corr = t->GetEntries(cut1)/static_cast<double>(t->GetEntries(cut2));
//	std::cout << "Correction factor: " << corr << std::endl;
//	f->Close();
//	return corr;
//}

bool MCJets::weightMC(jetType type, TH2D* effTrueHist, TH2D* effRecoHist) {
	if(!_inputsSet) {
		std::cout << "ERROR in MCJets::weightMC: input files not set yet" << std::endl;
		return false;
	}

	TString nameStr=typeName(type);

	TString histsFileName = gSaveDir + "/effHists_" + _name + "_" + nameStr + ".root"; 

	if(!_recreateInputs && !gSystem->AccessPathName(outputName(type))) {
		TFile* histsFile = TFile::Open(histsFileName);
		std::vector<TH2D*> hists {effTrueHist, effRecoHist};
		if(histsFile) {
			//we need to load both histograms
			//if one fails then we need to go through and recreate them
			bool good(true);

			if(effTrueHist) {
				TString hname = effTrueHist->GetName();
				TH2D* saved = static_cast<TH2D*>(histsFile->Get(hname));
				uint nbinsX = effTrueHist->GetNbinsX();
				uint nbinsY = effTrueHist->GetNbinsY();
				if(!saved ||
				   saved->GetNbinsX() != nbinsX || 
				   saved->GetNbinsY() != nbinsY) {
					good = false;
				} else {
					for(uint ibin=0; ibin<nbinsX; ++ibin) {
						for(uint jbin=0; jbin<nbinsY; ++jbin) {
							uint bin = effTrueHist->GetBin(ibin+1,jbin+1);
							effTrueHist->SetBinContent(bin,saved->GetBinContent(bin));
							effTrueHist->SetBinError(bin,saved->GetBinError(bin));
						}
					}
				}
			}
			if(effRecoHist) {
				TString hname = effRecoHist->GetName();
				TH2D* saved = static_cast<TH2D*>(histsFile->Get(hname));
				uint nbinsX = effRecoHist->GetNbinsX();
				uint nbinsY = effRecoHist->GetNbinsY();
				if(!saved ||
				   saved->GetNbinsX() != nbinsX || 
				   saved->GetNbinsY() != nbinsY) {
					good = false;
				} else {
					for(uint ibin=0; ibin<nbinsX; ++ibin) {
						for(uint jbin=0; jbin<nbinsY; ++jbin) {
							uint bin = effRecoHist->GetBin(ibin+1,jbin+1);
							effRecoHist->SetBinContent(bin,saved->GetBinContent(bin));
							effRecoHist->SetBinError(bin,saved->GetBinError(bin));
						}
					}
				}
			}
			histsFile->Close();
			if(good) {
				std::cout << "INFO in MCJets::weightMC: weighting already done" << std::endl;
				return true;
			}
		}

		//if we failed to get all histograms then reset any we did get
		if(effTrueHist) effTrueHist->Reset();
		if(effRecoHist) effRecoHist->Reset();
	}

	TH2D *effTrueNum(0), *effTrueDen(0);
	if(effTrueHist) {
		effTrueNum = cloneTH2D(TString(effTrueHist->GetName())+"_num",effTrueHist);
		effTrueDen = cloneTH2D(TString(effTrueHist->GetName())+"_den",effTrueHist);
		effTrueNum->Sumw2();
		effTrueDen->Sumw2();
	}
	TH2D *effRecoNum(0), *effRecoDen(0);
	if(effRecoHist) {
		effRecoNum = cloneTH2D(TString(effRecoHist->GetName())+"_num",effRecoHist);
		effRecoDen = cloneTH2D(TString(effRecoHist->GetName())+"_den",effRecoHist);
		effRecoNum->Sumw2();
		effRecoDen->Sumw2();
	}

	DatasetManager* dm = DatasetManager::getInstance();
//	TFile* fin(0);
	std::cout << "INFO : weighting " << nameStr << " MC sample" << std::endl;
	switch(type) {
		case jetAll0:
//			fin = TFile::Open(_lightInput);
			dm->setDataset(_lightInput);
			break;
		case jetAll4:
		case jetRecoD04:
		case jetRecoDp4:
		case jetRecoDs4:
		case jetRecoLc4:
		case jetRecoSV4:
		case jetTrueD04:
		case jetTrueDp4:
//			fin = TFile::Open(_charmInput);
			dm->setDataset(_charmInput);
			break;
		case jetAll5:
		case jetRecoD05:
		case jetRecoDp5:
		case jetRecoDs5:
		case jetRecoLc5:
		case jetRecoSV5:
		case jetTrueD05:
		case jetTrueDp5:
//			fin = TFile::Open(_beautyInput);
			dm->setDataset(_beautyInput);
			break;
		default:
			return false;
	}
	//if(!fin) return 0;
	//TTree* tin = dynamic_cast<TTree*>(fin->Get("T"));
	//if(!tin) return 0;

	double JetPT;
	double JetTruePT;
	double JetTrueEta;
	double JetTrueD0;
	double JetTrueDP;
	double JetTrueDSV;
	double JetTrueBSV;
	double JetTrueDPT;
	double JetTrueBPT;
	double NPV;
	double ptWeight(1.);

	double ZTRUEPZ, ZTRUEE, ZTRUEY;

	std::vector<double>* vSVMCor = new std::vector<double>();
	std::vector<double>* vSVN = new std::vector<double>();

	std::vector<double>* vD0M = new std::vector<double>();
	std::vector<double>* vD0PT = new std::vector<double>();
	std::vector<double>* vD0IPCHI2 = new std::vector<double>();
	std::vector<double>* vD0KWEIGHT = new std::vector<double>();
	std::vector<double>* vD0PIWEIGHT = new std::vector<double>();

	std::vector<double>* vDM = new std::vector<double>();
	std::vector<double>* vDPT = new std::vector<double>();
	std::vector<double>* vDIPCHI2 = new std::vector<double>();
	std::vector<double>* vDKWEIGHT = new std::vector<double>();
	std::vector<double>* vDPI1WEIGHT = new std::vector<double>();
	std::vector<double>* vDPI2WEIGHT = new std::vector<double>();

	dm->setBranchAddress("JetPT",       &JetPT);
	dm->setBranchAddress("JetTruePT",   &JetTruePT);
	dm->setBranchAddress("JetTrueEta",  &JetTrueEta);
	dm->setBranchAddress("JetTRUED0",   &JetTrueD0);
	dm->setBranchAddress("JetTRUEDP",   &JetTrueDP);
	dm->setBranchAddress("JetTRUEDSV",  &JetTrueDSV);
	dm->setBranchAddress("JetTRUEBSV",  &JetTrueBSV);
	dm->setBranchAddress("JetTRUEDPT",  &JetTrueDPT);
	dm->setBranchAddress("JetTRUEBPT",  &JetTrueBPT);

	dm->setBranchAddress("ZTRUEPZ",   &ZTRUEPZ);
	dm->setBranchAddress("ZTRUEE",   &ZTRUEE);

	dm->setBranchAddress("SVMCor",      &vSVMCor);
	dm->setBranchAddress("SVN",         &vSVN);

	dm->setBranchAddress("D0M",         &vD0M);
	dm->setBranchAddress("D0PT",        &vD0PT);
	dm->setBranchAddress("D0IPCHI2",    &vD0IPCHI2);
	dm->setBranchAddress("D0KWEIGHT",   &vD0KWEIGHT);
	dm->setBranchAddress("D0PIWEIGHT",  &vD0PIWEIGHT);

	dm->setBranchAddress("DM",         &vDM);
	dm->setBranchAddress("DPT",        &vDPT);
	dm->setBranchAddress("DIPCHI2",    &vDIPCHI2);
	dm->setBranchAddress("DKWEIGHT",   &vDKWEIGHT);
	dm->setBranchAddress("DPI1WEIGHT",  &vDPI1WEIGHT);
	dm->setBranchAddress("DPI2WEIGHT",  &vDPI2WEIGHT);

	dm->setBranchAddress("NPV",  &NPV);

	//weight MC for continuous true jet pT
	if(_useInputPtWeights) {
		std::cout << "INFO in MCJets:weightMC: using weight branch from input tuple as pT weights" << std::endl;
		dm->setBranchAddress("weight",  &ptWeight);
	} else if(_truePtWeights.find(type)==_truePtWeights.end()) {
		std::cout << "INFO in MCJets:weightMC: input true pT weights not set" << std::endl;
		std::cout << "                           using weight of 1" << std::endl;
		_truePtWeights[type] = new TH1D("jetTruePtWeights"+typeName(type),"",1,10e3,100e3);
		_truePtWeights[type]->SetBinContent(1,1.);
	}
	//end true(pT) weights setup

	TH1D* fPtDWeights = new TH1D("fPtDWeights","",20,0.,2.0);
	TH1D* ptRecWeights = new TH1D("ptRecWeights","",90,10.e3,100.e3);

	//TH1D* fPtDData = new TH1D("fPtDData"+typeName(type),"",20,0.,2.0);
	//TH1D* fPtDSim = new TH1D("fPtDSim"+typeName(type),"",20,0.,2.0);
	//TH1D* fPtDSB = new TH1D("fPtDSB"+typeName(type),"",20,0.,2.0);
		
	//for D0/Dp also weight to match pT(D0)/pT(jet) to data
	if(!_d0DataInput.IsNull() && (type==jetRecoD04 || type==jetRecoD05)) {
		TFile* fdata = TFile::Open(_d0DataInput);
		TTree* tdata = dynamic_cast<TTree*>(fdata->Get("T7"));

		TH1D* fPtDData = new TH1D("fPtDData"+typeName(type),"",20,0.,2.0);
		TH1D* fPtDSim = new TH1D("fPtDSim"+typeName(type),"",20,0.,2.0);
		TH1D* fPtDSB = new TH1D("fPtDSB"+typeName(type),"",20,0.,2.0);

		fPtDData->Sumw2();
		fPtDSim->Sumw2();
		fPtDSB->Sumw2();

		dm->draw("D0PT/JetPT","D0M>0 && JetTRUED0 && (D0KPNNK>0.2 || D0KPT>25000. || D0KP>500000.)",fPtDSim);//TODO no pion PID cut && D0PIPNNPI>0.1",fPtDSim);
		if(type==jetRecoD04) {
			tdata->Draw("D0PT/JetPT>>fPtDData"+typeName(type),"weight4*(TMath::Abs(D0M-1865)<30.)*(D0LogIPChi2<2.5)");
			tdata->Draw("D0PT/JetPT>>fPtDSB"+typeName(type),"weight4*(TMath::Abs(D0M-1865)>30. && TMath::Abs(D0M-1865)<60.)*(D0LogIPChi2<2.5)");
		} else if(type==jetRecoD05) {
			tdata->Draw("D0PT/JetPT>>fPtDData"+typeName(type),"weight5*(TMath::Abs(D0M-1865)<30.)*(D0LogIPChi2>2.5)");
			tdata->Draw("D0PT/JetPT>>fPtDSB"+typeName(type),"weight5*(TMath::Abs(D0M-1865)>30. && TMath::Abs(D0M-1865)<60.)*(D0LogIPChi2>2.5)");
		}
		std::cout << "!" << _d0DataInput << " " << tdata->GetEntries() << " " << fPtDData->Integral() << " " << fPtDSim->Integral() << " " << fPtDSB->Integral() << std::endl;//TODO

		fPtDData->Add(fPtDSB,-1.);

		fPtDData->Scale(1./fPtDData->Integral());
		fPtDSim->Scale(1./fPtDSim->Integral());

		fPtDWeights->Divide(fPtDData,fPtDSim);

		fPtDSim->SetLineColor(kRed);
		TCanvas c;
		fPtDData->Draw();
		fPtDSim->Draw("same");
		c.SaveAs(gSaveDir+"/D0fPtReweighting"+_name+"_"+nameStr+".pdf");

		fPtDWeights->Draw();
		c.SaveAs(gSaveDir+"/D0fPtWeights"+_name+"_"+nameStr+".pdf");

		fdata->Close();
	} else if(!_dpDataInput.IsNull() && (type==jetRecoDp4 || type==jetRecoDp5)) {
		TFile* fdata = TFile::Open(_dpDataInput);
		TTree* tdata = dynamic_cast<TTree*>(fdata->Get("T7"));

		TH1D* fPtDData = new TH1D("fPtDData"+typeName(type),"",20,0.,2.0);
		TH1D* fPtDSim = new TH1D("fPtDSim"+typeName(type),"",20,0.,2.0);
		TH1D* fPtDSB = new TH1D("fPtDSB"+typeName(type),"",20,0.,2.0);

		fPtDData->Sumw2();
		fPtDSim->Sumw2();
		fPtDSB->Sumw2();

		dm->draw("DPT/JetPT","DM>0 && JetTRUEDP && (DKPNNK>0.2 || DKPT>25000. || DKP>500000.)",fPtDSim);//TODO no pion PID cut && D0PIPNNPI>0.1",fPtDSim);
		if(type==jetRecoDp4) {
			tdata->Draw("DPT/JetPT>>fPtDData"+typeName(type),"weight4*(TMath::Abs(DM-1870)<30.)*(DLogIPChi2<2.5)");
			tdata->Draw("DPT/JetPT>>fPtDSB"+typeName(type),"weight4*(TMath::Abs(DM-1870)>30. && TMath::Abs(DM-1870)<60.)*(DLogIPChi2<2.5)");
		} else if(type==jetRecoDp5) {
			tdata->Draw("DPT/JetPT>>fPtDData"+typeName(type),"weight5*(TMath::Abs(DM-1870)<30.)*(DLogIPChi2>2.5)");
			tdata->Draw("DPT/JetPT>>fPtDSB"+typeName(type),"weight5*(TMath::Abs(DM-1870)>30. && TMath::Abs(DM-1870)<60.)*(DLogIPChi2>2.5)");
		}
		std::cout << "!" << _dpDataInput << " " << tdata->GetEntries() << " " << fPtDData->Integral() << " " << fPtDSim->Integral() << " " << fPtDSB->Integral() << std::endl;//TODO

		fPtDData->Add(fPtDSB,-1.);

		fPtDData->Scale(1./fPtDData->Integral());
		fPtDSim->Scale(1./fPtDSim->Integral());

		fPtDWeights->Divide(fPtDData,fPtDSim);

		fPtDSim->SetLineColor(kRed);
		TCanvas c;
		fPtDData->Draw();
		fPtDSim->Draw("same");
		c.SaveAs(gSaveDir+"/DpfPtReweighting"+_name+"_"+nameStr+".pdf");

		fPtDWeights->Draw();
		c.SaveAs(gSaveDir+"/DpfPtWeights"+_name+"_"+nameStr+".pdf");

		fdata->Close();
	}
	//end f(pT) weights setup

	//for SV also weight in pT_reco(jet) to data
	if(!_svDataInput.IsNull() && (type==jetRecoSV4 || type==jetRecoSV5)) {
		TFile* fdata = TFile::Open(_svDataInput);
		TTree* tdata = dynamic_cast<TTree*>(fdata->Get("T"));

		TH1D* ptRecData = new TH1D("ptRecData"+typeName(type),"",90,10.e3,100.e3);
		TH1D* ptRecSim  = new TH1D("ptRecSim"+typeName(type), "",90,10.e3,100.e3);

		ptRecData->Sumw2();
		ptRecSim->Sumw2();

		tdata->Draw("JetPT>>ptRecData"+typeName(type),"SVM[0]>0");
		//if(type==jetRecoSV4) {
		//	tin->Draw("JetPT>>ptRecSim","SVM[0]>0 && JetTRUEDSV");
		//} else if(type==jetRecoSV5) {
		//	tin->Draw("JetPT>>ptRecSim","SVM[0]>0 && JetTRUEBSV");
		//}
		int nEntries = dm->getEntries();
		boost::progress_display progress( nEntries );
		while(dm->getNext()) {
		//for(int iEntry=0; iEntry<nEntries; ++iEntry)
			++progress;
			//tin->GetEntry(iEntry);
			if(vSVMCor->size()<1) continue;
			if(type==jetRecoSV4 && !JetTrueDSV) continue;
			if(type==jetRecoSV5 && !JetTrueBSV) continue;
			double truePtWeight = _truePtWeights[type]->GetBinContent(_truePtWeights[type]->FindBin(JetTruePT));
			ptRecSim->Fill(JetPT,truePtWeight);
		}

		ptRecData->Scale(1./ptRecData->Integral());
		ptRecSim->Scale(1./ptRecSim->Integral());

		ptRecWeights->Divide(ptRecData,ptRecSim);

		ptRecSim->SetLineColor(kRed);
		TCanvas c;
		ptRecData->Draw();
		ptRecSim->Draw("same");
		c.SaveAs(gSaveDir+"/jetPtRecReweighting"+_name+"_"+nameStr+".pdf");

		ptRecWeights->Draw();
		c.SaveAs(gSaveDir+"/jetPtRecWeights"+_name+"_"+nameStr+".pdf");

		fdata->Close();
	}

	//now make the weighted output file
	TString fname=outputName(type);//_name+"MCjets";
	//if(dataIsResampledMC) fname+=file;
	//fname+="_"+nameStr+".root";

	TTree* tout = new TTree("T","");
	tout->SetDirectory(0);

	double D0M(0.), D0PT(0.), D0LogIPChi2(0.), DM(0.), DPT(0.), DLogIPChi2(0.), SVMCor(0.), SVN(0.), weight(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetTruePT",    &JetTruePT);
	tout->Branch("JetTrueEta",   &JetTrueEta);
	tout->Branch("ZTRUEPZ",      &ZTRUEPZ);
	tout->Branch("ZTRUEE",       &ZTRUEE);
	tout->Branch("ZTRUEY",       &ZTRUEY);
	tout->Branch("D0M",          &D0M);
	tout->Branch("D0PT",         &D0PT);
	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
	tout->Branch("DM",           &DM);
	tout->Branch("DPT",          &DPT);
	tout->Branch("DLogIPChi2",   &DLogIPChi2);
	tout->Branch("SVMCor",       &SVMCor);
	tout->Branch("SVN",          &SVN);
	tout->Branch("weight",       &weight);
	tout->Branch("NPV",          &NPV);

	int nEntries = dm->getEntries();
	boost::progress_display progress( nEntries );
	//for(int iEntry=0; iEntry<nEntries; ++iEntry)
	dm->rewind();
	while(dm->getNext()) {
		++progress;
		//tin->GetEntry(iEntry);

		if(vD0M ->size()>0) {
			D0M  = vD0M->at(0);
			D0PT = vD0PT->at(0);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(0));
		} else {
			D0M  = 0.;
			D0PT = 0.;
			D0LogIPChi2 = 0.;
		}
		if(vDM ->size()>0) {
			DM  = vDM->at(0);
			DPT = vDPT->at(0);
			DLogIPChi2 = TMath::Log(vDIPCHI2->at(0));
		} else {
			DM  = 0.;
			DPT = 0.;
			DLogIPChi2 = 0.;
		}
		if(vSVMCor ->size()>0) {
			SVMCor = vSVMCor->at(0);
			SVN = vSVN->at(0);
		} else {
			SVMCor = -1.;
			SVN = 0.;
		}

		if(_useInputPtWeights) weight = ptWeight;
		else weight = _truePtWeights[type]->GetBinContent(_truePtWeights[type]->FindBin(JetTruePT));

		//extra weights if we have them
		if(_extraPtWeights.find(type)!=_extraPtWeights.end()) {
			weight *= _extraPtWeights[type]->GetBinContent(_extraPtWeights[type]->FindBin(JetPT));
		}

		if(ZTRUEE>0.) ZTRUEY = 0.5*TMath::Log((ZTRUEE+ZTRUEPZ)/(ZTRUEE-ZTRUEPZ));
		else ZTRUEY = 3.;//default value for datasets without Z
		if(effTrueDen && JetPT>10e3) effTrueDen->Fill(JetTruePT,ZTRUEY,weight);
		if(effRecoDen && JetPT>10e3) effRecoDen->Fill(JetPT,ZTRUEY,weight);

		//only use matching jets
		switch(type) {
			case jetRecoD04:
				if(JetTrueDPT<5e3) continue;
				if(!JetTrueD0) continue;
				if(D0PT < _dptmin || (_dptmax!=-1 && D0PT > _dptmax) ) continue;
				break;
			case jetRecoD05:
				if(JetTrueBPT<5e3) continue;
				if(!JetTrueD0) continue;
				if(D0PT < _dptmin || (_dptmax!=-1 && D0PT > _dptmax) ) continue;
				break;
			case jetRecoDp4:
				if(JetTrueDPT<5e3) continue;
				if(!JetTrueDP) continue;
				if(DPT < _dptmin || (_dptmax!=-1 && DPT > _dptmax) ) continue;
				break;
			case jetRecoDp5:
				if(JetTrueBPT<5e3) continue;
				if(!JetTrueDP) continue;
				if(DPT < _dptmin || (_dptmax!=-1 && DPT > _dptmax) ) continue;
				break;
			case jetRecoSV4:
				if(JetTrueDPT<5e3) continue;
				if(!JetTrueDSV) continue;
				if(SVMCor<0) continue;
				break;
			case jetRecoSV5:
				if(JetTrueBPT<5e3) continue;
				if(!JetTrueBSV) continue;
				if(SVMCor<0) continue;
				break;
			case jetTrueD04:
				if(JetTrueDPT<5e3) continue;
				if(!JetTrueD0) continue;
				break;
			case jetTrueD05:
				if(JetTrueBPT<5e3) continue;
				if(!JetTrueD0) continue;
				break;
			case jetTrueDp4:
				if(JetTrueDPT<5e3) continue;
				if(!JetTrueDP) continue;
				break;
			case jetTrueDp5:
				if(JetTrueBPT<5e3) continue;
				if(!JetTrueDP) continue;
				break;
			case jetAll0:
			case jetAll4:
			case jetAll5:
				break;
			default:
				return 0;
		}

		if(effTrueNum) effTrueNum->Fill(JetTruePT,ZTRUEY,weight);
		if(effRecoNum) effRecoNum->Fill(JetPT,ZTRUEY,weight);

		if(type==jetRecoSV4||type==jetRecoSV5) {
			if(!_svDataInput.IsNull()) weight *= ptRecWeights->GetBinContent(ptRecWeights->FindBin(JetPT));
		}
		if(type==jetRecoD04||type==jetRecoD05) {
			if(!_d0DataInput.IsNull()) weight *= fPtDWeights->GetBinContent(fPtDWeights->FindBin(D0PT/JetPT));

			//apply PID weights
			weight *= vD0KWEIGHT->at(0);
			weight *= vD0PIWEIGHT->at(0);
		}
		if(type==jetRecoDp4||type==jetRecoDp5) {
			if(!_dpDataInput.IsNull()) weight *= fPtDWeights->GetBinContent(fPtDWeights->FindBin(DPT/JetPT));

			//apply PID weights
			weight *= vDKWEIGHT->at(0);
			weight *= vDPI1WEIGHT->at(0);
			weight *= vDPI2WEIGHT->at(0);
		}

		tout->Fill();
	}
	TFile* histsFile = TFile::Open(histsFileName,"RECREATE");
	if(effTrueHist) {
		effTrueHist->Divide(effTrueNum,effTrueDen);
		std::cout << "pT_true" << std::endl;
		for(int i=0; i<effTrueHist->GetNbinsX(); ++i) {
			std::cout << effTrueHist->GetXaxis()->GetBinLowEdge(i+1) << "-" << effTrueHist->GetXaxis()->GetBinUpEdge(i+1) << "\t" << effTrueHist->GetBinContent(i+1,1) << "\t" << effTrueNum->GetBinContent(i+1,1) << "\t" << effTrueDen->GetBinContent(i+1,1) << std::endl;
		}
		histsFile->cd();
		effTrueHist->Write();
	}
	if(effRecoHist) {
		effRecoHist->Divide(effRecoNum,effRecoDen);
		std::cout << "pT_reco" << std::endl;
		for(int i=0; i<effRecoHist->GetNbinsX(); ++i) {
			std::cout << effRecoHist->GetXaxis()->GetBinLowEdge(i+1) << "-" << effRecoHist->GetXaxis()->GetBinUpEdge(i+1) << "\t" << effRecoHist->GetBinContent(i+1,1) << "\t" << effRecoNum->GetBinContent(i+1,1) << "\t" << effRecoDen->GetBinContent(i+1,1) << std::endl;
		}
		histsFile->cd();
		effRecoHist->Write();
	}
	histsFile->Close();

	dm->reset();

	TFile* fDWeights = TFile::Open(gSaveDir + "/dPtWeights_" + _name + "_" + nameStr + ".root","RECREATE"); 
	fDWeights->cd();
	fPtDWeights->Write();
	fDWeights->Close();

	TFile* fout = TFile::Open(fname,"RECREATE");
	tout->SetDirectory(fout);
	tout->Write();
	tout->AutoSave();
	fout->Close();

	return true;
}

//TODO reweite to use DatasetManager
RooUnfoldResponse* MCJets::trainUnfold(TH1D* binning, jetType type, uint ybin) {
	if(!_ybins && ybin!=0) {
		std::cout << "INFO in MCJets::trainUnfolding : y binning not set" << std::endl;
		std::cout << "                                 will use full range of y" << std::endl;
		ybin=0;
	}
	if(_ybins && ybin>static_cast<uint>(_ybins->GetNbinsX())) {
		std::cout << "INFO in MCJets::trainUnfolding : y binning has fewer than " << ybin << " bins" << std::endl;
		std::cout << "                                 will use full range of y" << std::endl;
		ybin=0;
	}
	double ymin(0.), ymax(10.);
	if(ybin>0) {
		ymin = _ybins->GetBinLowEdge(ybin);
		ymax = _ybins->GetBinLowEdge(ybin+1);
	}
	TFile* f(0);
	TString nameStr = typeName(type);
	if(ybin!=0) {
		nameStr+="_";
		nameStr+=ybin;
	}
	std::cout << "INFO in MCJets::trainUnfolding : training unfolding for " << nameStr << std::endl;

	//TString fname=_name+"MCjets";
	//if(dataIsResampledMC) fname+=file;
	//fname+="_"+nameStr+".root";
	TString fname = outputName(type);

	f = TFile::Open(fname);
	if(!f) {
		std::cout << "INFO in MCJets::trainUnfolding: weighted input file not ready" << std::endl;
		std::cout << "                                will attempt to create file now" << std::endl;
		weightMC(type);
		f = TFile::Open(fname);
	}
	if(!f) return 0;
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return 0;

	RooUnfoldResponse* resp = new RooUnfoldResponse(binning,binning);

	int countSkip(0), countMiss(0), countFake(0), countReal(0);

	double JetPT;
	double JetScaledPT;
	double JetTruePT;
	double JetTrueEta;
	double ZY(3.);
	double SVMCor;
	double D0PT;
	double DPT;
	double weight(1.);

	t->SetBranchAddress("JetPT",       &JetPT);
	t->SetBranchAddress("JetTruePT",   &JetTruePT);
	t->SetBranchAddress("JetTrueEta",  &JetTrueEta);
	if(_ybins) {
		t->SetBranchAddress("ZTRUEY",  &ZY);
	}

	t->SetBranchAddress("SVMCor",      &SVMCor);
	t->SetBranchAddress("D0PT",        &D0PT);
	t->SetBranchAddress("DPT",         &DPT);
	t->SetBranchAddress("weight",      &weight);

	unsigned int nentries = t->GetEntries();
	boost::progress_display progress( nentries );
	for(unsigned int ientry=0; ientry<nentries; ++ientry) {
		++progress;
		t->GetEntry(ientry);
		if(JetTrueEta<2.2||JetTrueEta>4.2) continue;

		//apply jet energy scale + smearing
		JetScaledPT = gRandom->Gaus(_unfoldingScaleFactor*JetPT,
				           _unfoldingSmearFactor*_unfoldingScaleFactor*JetPT);

		//only use matching jets
		switch(type) {
			case jetRecoD04:
			case jetRecoD05:
				if(D0PT < _dptmin ||
				   (_dptmax!=-1 && D0PT > _dptmax) ) continue;
				break;
			case jetRecoDp4:
			case jetRecoDp5:
				if(DPT < _dptmin ||
				   (_dptmax!=-1 && DPT > _dptmax) ) continue;
				break;
			case jetRecoSV4:
				if(SVMCor<0) continue;
				break;
			case jetRecoSV5:
				if(SVMCor<0) continue;
				break;
			default:
				break;
		}
		if(ZY<ymin || ZY>ymax) continue;

		int truebin(0), recobin(0);
		int nbins = binning->GetNbinsX();

		truebin = binning->FindBin(JetTruePT);
		recobin = binning->FindBin(JetScaledPT);

		if((recobin<=0 || recobin>nbins) &&
		   (truebin<=0 || truebin>nbins) ) {
			++countSkip;
			continue;
		} else if(recobin<=0 || recobin>nbins) {
			++countMiss;
			resp->Miss(JetTruePT,weight);
//		} else if(truebin<=0 || truebin>nbins) {//TODO fakes break Ids method
//			++countFake;
//			resp->Fake(JetScaledPT,weight);
		} else {
			++countReal;
			resp->Fill(JetScaledPT,JetTruePT,weight);
		}
	}

	f->Close();

	std::cout << countReal << "\t" << countFake << "\t" << countMiss << "\t" << countSkip << std::endl;

	_unfoldings[type][ybin] = resp;

	return resp; //new RooUnfoldResponse(*resp);
}

TH1D* MCJets::unfold(TH1D* input, jetType type, uint ybin) {
	if(!input) return 0;

	TString nameStr=typeName(type);
	if(ybin!=0) {
		nameStr+="_";
		nameStr+=ybin;
	}
	std::cout << "INFO : unfolding " << nameStr << " results for file " << _name << std::endl;

	if(_dptmin!=5000. || _dptmax!=-1) {
		nameStr+="_D0";
		nameStr+=_dptmin;
		nameStr+="-";
		nameStr+=_dptmax;
	}

	//train unfolding on MC
	//use binning from data histogram
	RooUnfoldResponse* response(0);
	if(_unfoldings.find(type) == _unfoldings.end()) {
		_unfoldings[type] = std::map<uint,RooUnfoldResponse*>();
	}
	if(_unfoldings[type].find(ybin) == _unfoldings[type].end()) {
		response = trainUnfold(input,type,ybin);
	} else {
		response = _unfoldings[type][ybin];
	}

	if(!response) return 0;

	//plot response matrix
	TH1D* mcMeas  = dynamic_cast<TH1D*>(response->Hmeasured());
	TH1D* mcTrue  = dynamic_cast<TH1D*>(response->Htruth());
	TH2D* respMat = dynamic_cast<TH2D*>(response->Hresponse());

	setLHCbStyle();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 2;//5;
	const Int_t NCont = 500;
	Double_t stops[NRGBs]  = { 0.00, /*0.15, 0.40, 0.65,*/ 1.00};
	Double_t reds[NRGBs]   = { 1.00, /*1.00, 1.00, 1.00,*/ 0.00};
	Double_t greens[NRGBs] = { 1.00, /*0.95, 0.50, 0.00,*/ 0.00};
	Double_t blues[NRGBs]  = { 1.00, /*0.00, 0.00, 0.00,*/ 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);
	TCanvas c;
	respMat->GetXaxis()->SetTitle("#it{p}_{T}^{reco}(jet) [MeV/#it{c}]");
	respMat->GetYaxis()->SetTitle("#it{p}_{T}^{true}(jet) [MeV/#it{c}]");
	respMat->SetTitle("");
	respMat->Draw("colz");
	respMat->Draw("col");
	c.UseCurrentStyle();
	c.SetRightMargin(0.10);
	c.SetTopMargin(0.07);
	c.SaveAs(gSaveDir+"/unfoldingResponse_"+nameStr+".pdf");

	//apply response
	RooUnfold* unfolded(0);

	switch(_unfoldingMethod) {
		case unfoldInvert:
			unfolded = new RooUnfoldInvert(response, input);
			break;
		case unfoldBayes:
			if(_regPar>=0) unfolded = new RooUnfoldBayes(response, input, _regPar);
			else unfolded = new RooUnfoldBayes(response, input);
			break;
		case unfoldIDS:
			if(_regPar>=0) unfolded = new RooUnfoldIds(response, input, _regPar);
			else unfolded = new RooUnfoldIds(response, input);
			break;
		case unfoldSVD:
			if(_regPar>=0) unfolded = new RooUnfoldSvd(response, input, _regPar);
			else unfolded = new RooUnfoldSvd(response, input);
			break;
		case unfoldTUnfold:
			if(_regPar>=0) unfolded = new RooUnfoldTUnfold(response, input, _regPar);
			else unfolded = new RooUnfoldTUnfold(response, input);
			break;
		default:
			std::cout << "Unknown unfolding method " << _unfoldingMethod << std::endl;
			return 0;
	}
	//RooUnfoldDagostini unfolded(response, input, 1);
	TH1D* ret = dynamic_cast<TH1D*>(unfolded->Hreco());

	//TODO scaling messes up the output file
	//TODO//mcTrue->Scale(input->Integral()/mcMeas->Integral());
	//TODO//mcMeas->Scale(input->Integral()/mcMeas->Integral());

	//TODO//std::cout << nameStr <<" "<<input->GetMaximum()<<" "<<ret->GetMaximum()<<" "<<mcMeas->GetMaximum()<<" "<<mcTrue->GetMaximum() << std::endl;//TODO
	mcTrue->SetMarkerStyle(22);
	mcTrue->SetLineStyle(kDashed);
	input->SetLineColor(kRed);
	input->SetMarkerColor(kRed);
	ret->SetLineColor(kRed);
	ret->SetMarkerColor(kRed);
	ret->SetMarkerStyle(22);
	ret->SetLineStyle(kDashed);
	input->SetMaximum(1.1*TMath::Max(TMath::Max(input->GetMaximum(),ret->GetMaximum()),TMath::Max(mcMeas->GetMaximum()*input->Integral()/mcMeas->Integral(),mcTrue->GetMaximum()*input->Integral()/mcMeas->Integral())));
	input->SetMinimum(0.);
	input->GetXaxis()->SetTitle("#it{p}_{T}(jet) [MeV/#it{c}]");
	input->GetYaxis()->SetTitle("Candidates");
	input->SetTitle("");
	input->GetYaxis()->SetTitleOffset(1.1);
	input->Draw();
	mcMeas->DrawNormalized("same",input->Integral());
	mcTrue->DrawNormalized("same",input->Integral());
	input->Draw("same");
	ret->Draw("same");

	//TODO//mcMeas->SetMaximum(1.1*TMath::Max(TMath::Max(input->GetMaximum(),ret->GetMaximum()),TMath::Max(mcMeas->GetMaximum(),mcTrue->GetMaximum())));
	//TODO//mcMeas->SetMinimum(0.);
	//TODO//mcMeas->GetXaxis()->SetTitle("#it{p}_{T}(jet) [MeV/#it{c}]");
	//TODO//mcMeas->GetYaxis()->SetTitle("Candidates");
	//TODO//mcMeas->SetTitle("");
	//TODO//mcMeas->Draw();
	//TODO//mcTrue->Draw("same");
	//TODO//input->Draw("same");
	//TODO//ret->Draw("same");
	//TODO//c.UseCurrentStyle();
	//TODO//mcMeas->GetYaxis()->SetTitleOffset(1.1);
	//TODO//mcTrue->SetMarkerStyle(22);
	//TODO//mcTrue->SetLineStyle(kDashed);
	//TODO//input->SetLineColor(kRed);
	//TODO//input->SetMarkerColor(kRed);
	//TODO//ret->SetLineColor(kRed);
	//TODO//ret->SetMarkerColor(kRed);
	//TODO//ret->SetMarkerStyle(22);
	//TODO//ret->SetLineStyle(kDashed);
	c.SetRightMargin(0.10);
	c.SetTopMargin(0.07);
	c.SetLeftMargin(0.16);
	c.SaveAs(gSaveDir+"/unfoldingMC_"+nameStr+".pdf");

	//input->SetName("input");
	mcMeas->SetName("mcMeas");
	mcTrue->SetName("mcTrue");
	respMat->SetName("respMat");
	ret->SetName("ret");
	TFile* fsave = TFile::Open(gSaveDir+"/unfoldingMC_"+nameStr+".root","RECREATE");
	fsave->cd();
	input->Write();
	mcMeas->Write();
	mcTrue->Write();
	ret->Write();
	respMat->Write();

	return ret;
}

TH2D* MCJets::unfold(TH2D* input, jetType type, bool useYBins) {
	TH2D* unfolded = cloneTH2D(TString(input->GetName())+"_unfolded",input);
	if(useYBins) {
		TH1D* ybins = input->ProjectionY(TString(input->GetName())+"_ybins");
		setYBins(ybins);
	}
	for(int j=1; j<=input->GetNbinsY(); ++j) {
		TH1D* inputSlice = input->ProjectionX(TString(input->GetName())+"_slice",j,j);
		TH1D* unfldSlice(0);

		if(useYBins) unfldSlice = unfold(inputSlice,type,j);
		else unfldSlice = unfold(inputSlice,type);

		for(int i=1; i<=input->GetNbinsX(); ++i) {
			unfolded->SetBinContent(i,j,unfldSlice->GetBinContent(i));
			unfolded->SetBinError(i,j,unfldSlice->GetBinError(i));
		}
	}

	//run over y-integrated histogram for plots
	TH1D* inputInt = input->ProjectionX(TString(input->GetName())+"_slice");
	unfold(inputInt,type);

	return unfolded;
}

//TString MCJets::setupTrainTestSamples(int sample) {
//	if(!_inputsSet) {
//		std::cout << "ERROR in MCJets::setupTrainTestSamples: input files not set yet" << std::endl;
//		return "";
//	}
//	if(_mcSampled) {
//		std::cout << "WARNING in MCJets::setupTrainTestSamples: samples for MC study already constructed" << std::endl;
//		return _resampledName;
//	}
//
//	std::cout << "INFO in MCJets::setupTrainTestSamples: making samples for MC study " << sample << std::endl;
//	TRandom3 rand(1000+sample);
//	TString name = "_"; name+=sample;
//	TFile* fin4 = TFile::Open(_charmInput);
//	TTree* tin4 = dynamic_cast<TTree*>(fin4->Get("T"));
//
//	//TFile* fout40 = TFile::Open(gSaveDir+"/train_sim4"+name+".root","RECREATE");
//	//TTree* tout40 = tin4->CloneTree(0);
//
//	TFile* fout41 = TFile::Open(gSaveDir+"/test_sim4"+name+".root","RECREATE");
//	TTree* tout41 = tin4->CloneTree(0);
//
//	TFile* fin5 = TFile::Open(_beautyInput);
//	TTree* tin5 = dynamic_cast<TTree*>(fin5->Get("T"));
//
//	//TFile* fout50 = TFile::Open(gSaveDir+"/train_sim5"+name+".root","RECREATE");
//	//TTree* tout50 = tin5->CloneTree(0);
//
//	TFile* fout51 = TFile::Open(gSaveDir+"/test_sim5"+name+".root","RECREATE");
//	TTree* tout51 = tin5->CloneTree(0);
//
//	boost::progress_display progress( tin4->GetEntries()+tin5->GetEntries() );
//	for(int ientry=0; ientry<tin4->GetEntries(); ++ientry) {
//		++progress;
//		tin4->GetEntry(ientry);
//		if(rand.Rndm()<0.1) tout41->Fill();
//		//else tout41->Fill();
//	}
//
//	for(int ientry=0; ientry<tin5->GetEntries(); ++ientry) {
//		++progress;
//		tin5->GetEntry(ientry);
//		if(rand.Rndm()<0.1) tout51->Fill();
//		//else tout51->Fill();
//	}
//
//	//tout40->AutoSave();
//	//fout40->Close();
//	tout41->AutoSave();
//	fout41->Close();
//	//tout50->AutoSave();
//	//fout50->Close();
//	tout51->AutoSave();
//	fout51->Close();
//
//	//now merge the test samples
//	TChain* tout45 = new TChain("T");
//	tout45->Add(gSaveDir+"/test_sim4"+name+".root");
//	tout45->Add(gSaveDir+"/test_sim5"+name+".root");
//	tout45->Merge(gSaveDir+"/test_sim"+name+".root");
//
//	//update the input locations 
////	_charmInput = gSaveDir+"/train_sim4"+name+".root";
////	_beautyInput = gSaveDir+"/train_sim5"+name+".root";
//	_mcSampled=true;
//	_resampledName=gSaveDir+"/test_sim"+name+".root";
//
//	return _resampledName;
//}

//for backwards compatibility
//will put the histograms all into a map then pass them to the new version which is just going to unpack them again
bool MCJets::getTruth(TString file, TH1D* true4, TH1D* true5, TH1D* trueD04, TH1D* trueD05, TH1D* trueD0Sel4, TH1D* trueD0Sel5, TH1D* trueSV4, TH1D* trueSV5, bool useTruePT) {
	std::map<std::string, TH1D*> hists;
	if(true4) hists["true4"] = true4;
	if(true5) hists["true5"] = true5;
	if(trueD04) hists["trueD04"] = trueD04;
	if(trueD05) hists["trueD05"] = trueD05;
	if(trueD0Sel4) hists["trueD0Sel4"] = trueD0Sel4;
	if(trueD0Sel5) hists["trueD0Sel5"] = trueD0Sel5;
	if(trueSV4) hists["trueSV4"] = trueSV4;
	if(trueSV5) hists["trueSV5"] = trueSV5;

	return getTruth(file, hists, useTruePT);
}

bool MCJets::getTruth(TString file, std::map<std::string,TH1D*> hMap, bool useTruePT) {
	std::cout << "INFO in MCJets::getTruth: getting truth information" << std::endl;
	//if(!true4 || !true5 || !trueD04 || !trueD05 || !trueD0Sel4 || !trueD0Sel5 || !trueSV4 || !trueSV5) return false;
	//
	TH1D *true4(0), *true5(0), *trueD04(0), *trueD05(0), *trueD0Sel4(0), *trueD0Sel5(0), *trueDp4(0), *trueDp5(0), *trueDpSel4(0), *trueDpSel5(0), *trueSV4(0), *trueSV5(0);

	if(hMap.find("true4") != hMap.end()) true4 = hMap["true4"];
	if(hMap.find("true5") != hMap.end()) true5 = hMap["true5"];
	if(hMap.find("trueD04") != hMap.end()) trueD04 = hMap["trueD04"];
	if(hMap.find("trueD05") != hMap.end()) trueD05 = hMap["trueD05"];
	if(hMap.find("trueD0Sel4") != hMap.end()) trueD0Sel4 = hMap["trueD0Sel4"];
	if(hMap.find("trueD0Sel5") != hMap.end()) trueD0Sel5 = hMap["trueD0Sel5"];
	if(hMap.find("trueDp4") != hMap.end()) trueDp4 = hMap["trueDp4"];
	if(hMap.find("trueDp5") != hMap.end()) trueDp5 = hMap["trueDp5"];
	if(hMap.find("trueDpSel4") != hMap.end()) trueDpSel4 = hMap["trueDpSel4"];
	if(hMap.find("trueDpSel5") != hMap.end()) trueDpSel5 = hMap["trueDpSel5"];
	if(hMap.find("trueSV4") != hMap.end()) trueSV4 = hMap["trueSV4"];
	if(hMap.find("trueSV5") != hMap.end()) trueSV5 = hMap["trueSV5"];

	//check if we already calculated these
	TString histsFileName = gSaveDir+"/truthHists_"+_name+"_"+file;
	if(useTruePT) histsFileName+="_truePt";
	histsFileName+=".root";
	
	TFile* histsFile = TFile::Open(histsFileName);
	//container of all histograms to iterate
	std::vector<TH1D*> hists {true4,true5,trueD04,trueD05,trueD0Sel4,trueD0Sel5,trueDp4,trueDp5,trueDpSel4,trueDpSel5,trueSV4,trueSV5};
	if(histsFile) {
		//we need to load every histogram
		//if one fails then we need to go through and recreate them
		bool good(true);
		for(auto hist: hists) {
			if(!hist) continue;
			TString name = hist->GetName();
			TH1D* saved = static_cast<TH1D*>(histsFile->Get(name));
			if(!saved) {
				good = false;
				break;
			}

			uint nbins = hist->GetNbinsX();
			if(saved->GetNbinsX()!=nbins) {
				good = false;
				break;
			}
			for(uint ibin=0; ibin<nbins; ++ibin) {
				hist->SetBinContent(ibin+1,saved->GetBinContent(ibin+1));
				hist->SetBinError(ibin+1,saved->GetBinError(ibin+1));
			}

		}
		histsFile->Close();
		if(good) {
			std::cout << "INFO in MCJets::getTruth: truth histograms already created" << std::endl;
			return true;
		}
	}

	for(auto hist: hists) {
		if(!hist) continue;
		hist->Reset();
	}

	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset(file);
	//TFile* f = TFile::Open(file);
	//if(!f) return false;
	//TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	//if(!t) return false;

	double JetPT;
	double JetTruePT;
	double JetTrueD0;
	double JetTrueDP;
	double JetTrueDSV;
	double JetTrueBSV;
	double JetTRUEc;
	double JetTRUEb;
	double JetTrueDPT;
	double JetTrueBPT;

	std::vector<double>* TRUEDID = new std::vector<double>();
	std::vector<double>* TRUEDFROMB = new std::vector<double>();

	std::vector<double>* SVN = new std::vector<double>();
	std::vector<double>* D0M = new std::vector<double>();
	std::vector<double> *D0IPCHI2 = new std::vector<double>();
	std::vector<double> *D0PT = new std::vector<double>();
	std::vector<double> *D0PX = new std::vector<double>();
	std::vector<double> *D0PY = new std::vector<double>();
	std::vector<double> *D0PZ = new std::vector<double>();
	std::vector<double> *D0E = new std::vector<double>();
	std::vector<double> *D0X = new std::vector<double>();
	std::vector<double> *D0Y = new std::vector<double>();
	std::vector<double> *D0Z = new std::vector<double>();
	std::vector<double> *D0KP = new std::vector<double>();
	std::vector<double> *D0KPT = new std::vector<double>();
	std::vector<double> *D0KPX = new std::vector<double>();
	std::vector<double> *D0KPY = new std::vector<double>();
	std::vector<double> *D0KPZ = new std::vector<double>();
	std::vector<double> *D0PIP = new std::vector<double>();
	std::vector<double> *D0PIPT = new std::vector<double>();
	std::vector<double> *D0PIPX = new std::vector<double>();
	std::vector<double> *D0PIPY = new std::vector<double>();
	std::vector<double> *D0PIPZ = new std::vector<double>();
	std::vector<double> *D0KPNNK = new std::vector<double>();
	//std::vector<double> *D0PIPNNPI = new std::vector<double>();
	std::vector<double>* D0TRUEIDX = new std::vector<double>();
	std::vector<double>* DM = new std::vector<double>();
	std::vector<double> *DIPCHI2 = new std::vector<double>();
	std::vector<double> *DPT = new std::vector<double>();
	std::vector<double> *DPX = new std::vector<double>();
	std::vector<double> *DPY = new std::vector<double>();
	std::vector<double> *DPZ = new std::vector<double>();
	std::vector<double> *DE = new std::vector<double>();
	std::vector<double> *DX = new std::vector<double>();
	std::vector<double> *DY = new std::vector<double>();
	std::vector<double> *DZ = new std::vector<double>();
	std::vector<double> *DKP = new std::vector<double>();
	std::vector<double> *DKPT = new std::vector<double>();
	std::vector<double> *DKPX = new std::vector<double>();
	std::vector<double> *DKPY = new std::vector<double>();
	std::vector<double> *DKPZ = new std::vector<double>();
	std::vector<double> *DPI1P = new std::vector<double>();
	std::vector<double> *DPI1PT = new std::vector<double>();
	std::vector<double> *DPI1PX = new std::vector<double>();
	std::vector<double> *DPI1PY = new std::vector<double>();
	std::vector<double> *DPI1PZ = new std::vector<double>();
	std::vector<double> *DPI2P = new std::vector<double>();
	std::vector<double> *DPI2PT = new std::vector<double>();
	std::vector<double> *DPI2PX = new std::vector<double>();
	std::vector<double> *DPI2PY = new std::vector<double>();
	std::vector<double> *DPI2PZ = new std::vector<double>();
	std::vector<double> *DKPNNK = new std::vector<double>();
	std::vector<double>* DTRUEIDX = new std::vector<double>();

	dm->setBranchAddress("JetPT",      &JetPT);
	dm->setBranchAddress("JetTruePT",  &JetTruePT);
	dm->setBranchAddress("JetTRUED0",  &JetTrueD0);
	dm->setBranchAddress("JetTRUEDP",  &JetTrueDP);
	dm->setBranchAddress("JetTRUEDSV", &JetTrueDSV);
	dm->setBranchAddress("JetTRUEBSV", &JetTrueBSV);
	dm->setBranchAddress("JetTRUEc",   &JetTRUEc);
	dm->setBranchAddress("JetTRUEb",   &JetTRUEb);
	dm->setBranchAddress("JetTRUEDPT", &JetTrueDPT);
	dm->setBranchAddress("JetTRUEBPT", &JetTrueBPT);

	dm->setBranchAddress("TRUEDID",    &TRUEDID);
	dm->setBranchAddress("TRUEDFROMB", &TRUEDFROMB);

	dm->setBranchAddress("SVN",        &SVN);
	dm->setBranchAddress("D0M",        &D0M);
	dm->setBranchAddress("D0IPCHI2",   &D0IPCHI2);
	dm->setBranchAddress("D0PT",       &D0PT);
	dm->setBranchAddress("D0PX",       &D0PX);
	dm->setBranchAddress("D0PY",       &D0PY);
	dm->setBranchAddress("D0PZ",       &D0PZ);
	dm->setBranchAddress("D0E",        &D0E);
	dm->setBranchAddress("D0X",        &D0X);
	dm->setBranchAddress("D0Y",        &D0Y);
	dm->setBranchAddress("D0Z",        &D0Z);
	dm->setBranchAddress("D0KP",       &D0KP);
	dm->setBranchAddress("D0KPT",      &D0KPT);
	dm->setBranchAddress("D0KPX",      &D0KPX);
	dm->setBranchAddress("D0KPY",      &D0KPY);
	dm->setBranchAddress("D0KPZ",      &D0KPZ);
	dm->setBranchAddress("D0PIP",      &D0PIP);
	dm->setBranchAddress("D0PIPT",     &D0PIPT);
	dm->setBranchAddress("D0PIPX",     &D0PIPX);
	dm->setBranchAddress("D0PIPY",     &D0PIPY);
	dm->setBranchAddress("D0PIPZ",     &D0PIPZ);
	dm->setBranchAddress("D0KPNNK",    &D0KPNNK);
	//dm->setBranchAddress("D0PIPNNPI",  &D0PIPNNPI);
	dm->setBranchAddress("DTRUEIDX",  &DTRUEIDX);
	dm->setBranchAddress("DM",        &DM);
	dm->setBranchAddress("DIPCHI2",   &DIPCHI2);
	dm->setBranchAddress("DPT",       &DPT);
	dm->setBranchAddress("DPX",       &DPX);
	dm->setBranchAddress("DPY",       &DPY);
	dm->setBranchAddress("DPZ",       &DPZ);
	dm->setBranchAddress("DE",        &DE);
	dm->setBranchAddress("DX",        &DX);
	dm->setBranchAddress("DY",        &DY);
	dm->setBranchAddress("DZ",        &DZ);
	dm->setBranchAddress("DKP",       &DKP);
	dm->setBranchAddress("DKPT",      &DKPT);
	dm->setBranchAddress("DKPX",      &DKPX);
	dm->setBranchAddress("DKPY",      &DKPY);
	dm->setBranchAddress("DKPZ",      &DKPZ);
	dm->setBranchAddress("DPI1P",      &DPI1P);
	dm->setBranchAddress("DPI1PT",     &DPI1PT);
	dm->setBranchAddress("DPI1PX",     &DPI1PX);
	dm->setBranchAddress("DPI1PY",     &DPI1PY);
	dm->setBranchAddress("DPI1PZ",     &DPI1PZ);
	dm->setBranchAddress("DPI2P",      &DPI2P);
	dm->setBranchAddress("DPI2PT",     &DPI2PT);
	dm->setBranchAddress("DPI2PX",     &DPI2PX);
	dm->setBranchAddress("DPI2PY",     &DPI2PY);
	dm->setBranchAddress("DPI2PZ",     &DPI2PZ);
	dm->setBranchAddress("DKPNNK",    &DKPNNK);
	dm->setBranchAddress("DTRUEIDX",  &DTRUEIDX);

	boost::progress_display progress( dm->getEntries() );
	//for(int i=0; i<t->GetEntries(); ++i)
	while(dm->getNext()) {
		++progress;
		//t->GetEntry(i);
		if(JetPT<10e3) continue;//only use jets that have been reconstructed
		if(useTruePT) JetPT = JetTruePT;

		double weight(1.);
		//if(JetTruePT>50000.) {
		//	weight=0.007;
		//} else if(JetTruePT>20000.) {
		//	weight=0.10;
		//} else if(JetTruePT>15000.) {
		//	weight=0.25;
		//}

		//b-jets
		if(JetTRUEb && JetTrueBPT>5e3) {
			if(true5) true5->Fill(JetPT,weight);
			//b->D0
			if(JetTrueD0) {
				for(unsigned int id=0; id<TRUEDID->size(); ++id) {
					if(TMath::Abs(TRUEDID->at(id))!=421.) continue;
					if(TRUEDFROMB->at(id)) {
						if(trueD05) trueD05->Fill(JetPT,weight);
						//TODO//if(D0M->size()>0) trueD0Sel5->Fill(JetPT,weight);
						//TODO//if(D0M->size()>0 && D0TRUEIDX->at(0) == id) trueD0Sel5->Fill(JetPT,weight);
						for(unsigned int irec=0; irec<D0TRUEIDX->size(); ++irec) {
							TVector3 D0P (D0PX->at(irec)  ,D0PY->at(irec)  ,D0PZ->at(irec));
							TVector3 D0P0(D0KPX->at(irec) ,D0KPY->at(irec) ,D0KPZ->at(irec));
							TVector3 D0P1(D0PIPX->at(irec),D0PIPY->at(irec),D0PIPZ->at(irec));

							if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
							if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
							if(!(D0P.Pt()>5e3)) continue;

							if(D0KPNNK->at(irec)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT

							if(trueD0Sel5 && D0TRUEIDX->at(irec) == id) trueD0Sel5->Fill(JetPT,weight);
						}
						//Note: the actual requirement in data is somewhere in between these last two options
						// We only keep one D0 candidate but we apply some cuts before making that choice so 
						// it's not always at index 0
						break;
					}
				}
			}
			//b->D+
			if(JetTrueDP) {
				for(unsigned int id=0; id<TRUEDID->size(); ++id) {
					if(TMath::Abs(TRUEDID->at(id))!=411.) continue;
					if(TRUEDFROMB->at(id)) {
						if(trueDp5) trueDp5->Fill(JetPT,weight);
						for(unsigned int irec=0; irec<DTRUEIDX->size(); ++irec) {
							TVector3 DP (DPX->at(irec)   ,DPY->at(irec)   ,DPZ->at(irec));
							TVector3 DP0(DKPX->at(irec)  ,DKPY->at(irec)  ,DKPZ->at(irec));
							TVector3 DP1(DPI1PX->at(irec),DPI1PY->at(irec),DPI1PZ->at(irec));
							TVector3 DP2(DPI2PX->at(irec),DPI2PY->at(irec),DPI2PZ->at(irec));

							if(!(DP0.Eta()>2. && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
							if(!(DP1.Eta()>2. && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
							if(!(DP2.Eta()>2. && DP2.Eta()<4.5 && DP2.Pt()>500. && DP2.Mag()>5000.)) continue;
							if(!(DP.Pt()>5e3)) continue;

							if(DKPNNK->at(irec)<0.2 && DP0.Pt()<25000. && DP0.Mag()<500000.) continue; //kaon PID turned off for high p or pT

							if(trueDpSel5 && DTRUEIDX->at(irec) == id) trueDpSel5->Fill(JetPT,weight);
						}
						//Note: the actual requirement in data is somewhere in between these last two options
						// We only keep one D+ candidate but we apply some cuts before making that choice so 
						// it's not always at index 0
						break;
					}
				}
			}
			//b->SV
			if(true) {
			//if(JetTrueBSV) { //Note: turned off the requirement for the SV to be "real"
				if(trueSV5 && SVN->size()>0) trueSV5->Fill(JetPT,weight);
			}
		//c-jets
		} else if(JetTRUEc && JetTrueDPT>5e3) {
			if(true4) true4->Fill(JetPT,weight);
			//c->D0
			if(JetTrueD0) {
				for(unsigned int id=0; id<TRUEDID->size(); ++id) {
					if(TMath::Abs(TRUEDID->at(id))!=421.) continue;
					if(!TRUEDFROMB->at(id)) {
						if(trueD04) trueD04->Fill(JetPT,weight);
						//TODO//if(D0M->size()>0) trueD0Sel4->Fill(JetPT,weight);
						//TODO//if(D0M->size()>0 && D0TRUEIDX->at(0) == id) trueD0Sel4->Fill(JetPT,weight);
						for(unsigned int irec=0; irec<D0TRUEIDX->size(); ++irec) {
							TVector3 D0P (D0PX->at(irec)  ,D0PY->at(irec)  ,D0PZ->at(irec));
							TVector3 D0P0(D0KPX->at(irec) ,D0KPY->at(irec) ,D0KPZ->at(irec));
							TVector3 D0P1(D0PIPX->at(irec),D0PIPY->at(irec),D0PIPZ->at(irec));

							if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
							if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
							if(!(D0P.Pt()>5e3)) continue;

							if(D0KPNNK->at(irec)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT

							if(trueD0Sel4 && D0TRUEIDX->at(irec) == id) trueD0Sel4->Fill(JetPT,weight);
						}
						//Note: the actual requirement in data is somewhere in between these last two options
						// We only keep one D0 candidate but we apply some cuts before making that choice so 
						// it's not always at index 0
						break;
					}
				}
			}
			//c->D+
			if(JetTrueDP) {
				for(unsigned int id=0; id<TRUEDID->size(); ++id) {
					if(TMath::Abs(TRUEDID->at(id))!=411.) continue;
					if(!TRUEDFROMB->at(id)) {
						if(trueDp4) trueDp4->Fill(JetPT,weight);
						for(unsigned int irec=0; irec<DTRUEIDX->size(); ++irec) {
							TVector3 DP (DPX->at(irec)   ,DPY->at(irec)   ,DPZ->at(irec));
							TVector3 DP0(DKPX->at(irec)  ,DKPY->at(irec)  ,DKPZ->at(irec));
							TVector3 DP1(DPI1PX->at(irec),DPI1PY->at(irec),DPI1PZ->at(irec));
							TVector3 DP2(DPI2PX->at(irec),DPI2PY->at(irec),DPI2PZ->at(irec));

							if(!(DP0.Eta()>2. && DP0.Eta()<4.5 && DP0.Pt()>500. && DP0.Mag()>5000.)) continue;
							if(!(DP1.Eta()>2. && DP1.Eta()<4.5 && DP1.Pt()>500. && DP1.Mag()>5000.)) continue;
							if(!(DP2.Eta()>2. && DP2.Eta()<4.5 && DP2.Pt()>500. && DP2.Mag()>5000.)) continue;
							if(!(DP.Pt()>5e3)) continue;

							if(DKPNNK->at(irec)<0.2 && DP0.Pt()<25000. && DP0.Mag()<500000.) continue; //kaon PID turned off for high p or pT

							if(trueDpSel4 && DTRUEIDX->at(irec) == id) trueDpSel4->Fill(JetPT,weight);
						}
						//Note: the actual requirement in data is somewhere in between these last two options
						// We only keep one D+ candidate but we apply some cuts before making that choice so 
						// it's not always at index 0
						break;
					}
				}
			}
			//c->SV
			if(true) {
			//if(JetTrueDSV) { //Note: turned off the requirement for the SV to be "real"
				if(trueSV4 && SVN->size()>0) trueSV4->Fill(JetPT,weight);
			}
		}
	}

	histsFile = TFile::Open(histsFileName, "RECREATE");
	histsFile->cd();
	if(true4) true4->Write();
	if(true5) true5->Write();
	if(trueD04) trueD04->Write();
	if(trueD05) trueD05->Write();
	if(trueD0Sel4) trueD0Sel4->Write();
	if(trueD0Sel5) trueD0Sel5->Write();
	if(trueDp4) trueDp4->Write();
	if(trueDp5) trueDp5->Write();
	if(trueDpSel4) trueDpSel4->Write();
	if(trueDpSel5) trueDpSel5->Write();
	if(trueSV4) trueSV4->Write();
	if(trueSV5) trueSV5->Write();
	histsFile->Close();

	dm->reset();

	return true;
}

//2D variant of getTruth for Z+j fits
bool MCJets::getTruth(TString file, TH2D* true0, TH2D* true4, TH2D* true5, TH2D* trueSV4, TH2D* trueSV5, bool useTruePT) {
	std::cout << "INFO in MCJets::getTruth: getting truth information" << std::endl;
	if(!true0 || !true4 || !true5 || !trueSV4 || !trueSV5) return false;

	//check if we already calculated these
	TString histsFileName = gSaveDir+"/truthHists2D_"+_name+"_"+file;
	if(useTruePT) histsFileName+="_truePt";
	histsFileName+=".root";
	
	TFile* histsFile = TFile::Open(histsFileName);
	//container of all histograms to iterate
	std::vector<TH2D*> hists {true0,true4,true5,trueSV4,trueSV5};
	if(histsFile) {
		//we need to load every histogram
		//if one fails then we need to go through and recreate them
		bool good(true);
		for(auto hist: hists) {
			TString name = hist->GetName();
			TH2D* saved = static_cast<TH2D*>(histsFile->Get(name));
			if(!saved) {
				good = false;
				break;
			}

			uint nbins = hist->GetNbinsX();
			if(saved->GetNbinsX()!=nbins) {
				good = false;
				break;
			}
			for(uint ibin=0; ibin<nbins; ++ibin) {
				hist->SetBinContent(ibin+1,saved->GetBinContent(ibin+1));
				hist->SetBinError(ibin+1,saved->GetBinError(ibin+1));
			}

		}
		histsFile->Close();
		if(good) {
			std::cout << "INFO in MCJets::getTruth: truth histograms already created" << std::endl;
			return true;
		}
	}

	for(auto hist: hists) {
		hist->Reset();
	}

	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset(file);

	double JetPT;
	double JetTruePT;
	double JetTrueD0;
	double JetTrueDSV;
	double JetTrueBSV;
	double JetTRUEc;
	double JetTRUEb;
	double JetTrueDPT;
	double JetTrueBPT;

	double ZTRUEPZ;
	double ZTRUEE;

	std::vector<double>* TRUEDID = new std::vector<double>();
	std::vector<double>* TRUEDFROMB = new std::vector<double>();

	std::vector<double>* SVN = new std::vector<double>();
	std::vector<double>* D0M = new std::vector<double>();
	std::vector<double> *D0IPCHI2 = new std::vector<double>();
	std::vector<double> *D0PT = new std::vector<double>();
	std::vector<double> *D0PX = new std::vector<double>();
	std::vector<double> *D0PY = new std::vector<double>();
	std::vector<double> *D0PZ = new std::vector<double>();
	std::vector<double> *D0E = new std::vector<double>();
	std::vector<double> *D0X = new std::vector<double>();
	std::vector<double> *D0Y = new std::vector<double>();
	std::vector<double> *D0Z = new std::vector<double>();
	std::vector<double> *D0KP = new std::vector<double>();
	std::vector<double> *D0KPT = new std::vector<double>();
	std::vector<double> *D0KPX = new std::vector<double>();
	std::vector<double> *D0KPY = new std::vector<double>();
	std::vector<double> *D0KPZ = new std::vector<double>();
	std::vector<double> *D0PIP = new std::vector<double>();
	std::vector<double> *D0PIPT = new std::vector<double>();
	std::vector<double> *D0PIPX = new std::vector<double>();
	std::vector<double> *D0PIPY = new std::vector<double>();
	std::vector<double> *D0PIPZ = new std::vector<double>();
	std::vector<double> *D0KPNNK = new std::vector<double>();
	std::vector<double> *D0PIPNNPI = new std::vector<double>();
	std::vector<double>* D0TRUEIDX = new std::vector<double>();

	dm->setBranchAddress("JetPT",      &JetPT);
	dm->setBranchAddress("JetTruePT",  &JetTruePT);
	dm->setBranchAddress("JetTRUED0",  &JetTrueD0);
	dm->setBranchAddress("JetTRUEDSV", &JetTrueDSV);
	dm->setBranchAddress("JetTRUEBSV", &JetTrueBSV);
	dm->setBranchAddress("JetTRUEc",   &JetTRUEc);
	dm->setBranchAddress("JetTRUEb",   &JetTRUEb);
	dm->setBranchAddress("JetTRUEDPT", &JetTrueDPT);
	dm->setBranchAddress("JetTRUEBPT", &JetTrueBPT);

	dm->setBranchAddress("ZTRUEPZ",    &ZTRUEPZ);
	dm->setBranchAddress("ZTRUEE",     &ZTRUEE);

	dm->setBranchAddress("TRUEDID",    &TRUEDID);
	dm->setBranchAddress("TRUEDFROMB", &TRUEDFROMB);

	dm->setBranchAddress("SVN",        &SVN);
	dm->setBranchAddress("D0M",        &D0M);
	dm->setBranchAddress("D0IPCHI2",   &D0IPCHI2);
	dm->setBranchAddress("D0PT",       &D0PT);
	dm->setBranchAddress("D0PX",       &D0PX);
	dm->setBranchAddress("D0PY",       &D0PY);
	dm->setBranchAddress("D0PZ",       &D0PZ);
	dm->setBranchAddress("D0E",        &D0E);
	dm->setBranchAddress("D0X",        &D0X);
	dm->setBranchAddress("D0Y",        &D0Y);
	dm->setBranchAddress("D0Z",        &D0Z);
	dm->setBranchAddress("D0KP",       &D0KP);
	dm->setBranchAddress("D0KPT",      &D0KPT);
	dm->setBranchAddress("D0KPX",      &D0KPX);
	dm->setBranchAddress("D0KPY",      &D0KPY);
	dm->setBranchAddress("D0KPZ",      &D0KPZ);
	dm->setBranchAddress("D0PIP",      &D0PIP);
	dm->setBranchAddress("D0PIPT",     &D0PIPT);
	dm->setBranchAddress("D0PIPX",     &D0PIPX);
	dm->setBranchAddress("D0PIPY",     &D0PIPY);
	dm->setBranchAddress("D0PIPZ",     &D0PIPZ);
	dm->setBranchAddress("D0KPNNK",    &D0KPNNK);
	dm->setBranchAddress("D0PIPNNPI",  &D0PIPNNPI);
	dm->setBranchAddress("D0TRUEIDX",  &D0TRUEIDX);

	boost::progress_display progress( dm->getEntries() );
	//for(int i=0; i<t->GetEntries(); ++i)
	while(dm->getNext()) {
		++progress;
		//t->GetEntry(i);
		//TODO off for now//if(JetPT<10e3) continue;//only use jets that have been reconstructed
		if(useTruePT) JetPT = JetTruePT;

		double ZEta = 0.5*TMath::Log((ZTRUEE+ZTRUEPZ)/(ZTRUEE-ZTRUEPZ));

		double weight(1.);
		//if(JetTruePT>50000.) {
		//	weight=0.007;
		//} else if(JetTruePT>20000.) {
		//	weight=0.10;
		//} else if(JetTruePT>15000.) {
		//	weight=0.25;
		//}

		true0->Fill(JetPT,ZEta,weight);
		//b-jets
		if(JetTRUEb && JetTrueBPT>5e3) {
			true5->Fill(JetPT,ZEta,weight);
			//b->SV
			if(true) {
			//if(JetTrueBSV) { //Note: turned off the requirement for the SV to be "real"
				if(SVN->size()>0) trueSV5->Fill(JetPT,ZEta,weight);
			}
		//c-jets
		} else if(JetTRUEc && JetTrueDPT>5e3) {
			true4->Fill(JetPT,ZEta,weight);
			//c->SV
			if(true) {
			//if(JetTrueDSV) { //Note: turned off the requirement for the SV to be "real"
				if(SVN->size()>0) trueSV4->Fill(JetPT,ZEta,weight);
			}
		}
	}

	histsFile = TFile::Open(histsFileName, "RECREATE");
	histsFile->cd();
	true0->Write();
	true4->Write();
	true5->Write();
	trueSV4->Write();
	trueSV5->Write();
	histsFile->Close();

	dm->reset();

	return true;
}

const TString MCJets::outputName(jetType type) {
	return gSaveDir+"/"+_name+"MCjets_"+typeName(type)+".root";
}

const TString MCJets::typeName(jetType type) {
	switch(type) {
		case jetAll0:
			return "q2all";
		case jetAll4:
			return "c2all";
		case jetAll5:
			return "b2all";
		case jetRecoD04:
			return "c2D0";
		case jetRecoDp4:
			return "c2Dp";
		case jetRecoSV4:
			return "c2SV";
		case jetTrueD04:
			return "c2allD0";
		case jetTrueDp4:
			return "c2allDp";
		case jetRecoD05:
			return "b2D0";
		case jetRecoDp5:
			return "b2Dp";
		case jetRecoSV5:
			return "b2SV";
		case jetTrueD05:
			return "b2allD0";
		case jetTrueDp5:
			return "b2allDp";
		default:
			return "";
	}
}
