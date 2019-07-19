#include "MCJets.h"

#include <iostream>
#include <vector>

#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TTree.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"

#include "cloneHists.h"
#include "outputFunctions.h"

void MCJets::setInputs(TString light, TString charm, TString beauty, TString svData, TString d0Data) {
	_lightInput  = light;
	_charmInput  = charm;
	_beautyInput = beauty;
	_svDataInput = svData;
	_d0DataInput = d0Data;

	_inputsSet = true;
}

void MCJets::setInputTruePtWeights(jetType type, TH1D* truePtWeights) {
	_truePtWeights[type] = truePtWeights;
}

void MCJets::setD0PtRange(double minpt, double maxpt) {
	_d0ptmin = minpt;
	_d0ptmax = maxpt;
}

double MCJets::getPtCorrFactor(jetType type, double ptMin, double ptMax) {
	TFile* f(0);
	TTree* t(0);
	TString cut1;
	TString cut2;
	double corr(1.);

	cut1 = "JetPT>"; cut1+=ptMin; cut1+=" && JetPT<"; cut1+=ptMax;
	cut2 = "JetPT>"; cut2+=ptMin; cut2+=" && JetPT<"; cut2+=ptMax;

	switch(type) {
		case jetAll0:
			std::cout << "getting pT correction factor for all jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetAll4:
			std::cout << "getting pT correction factor for all c-jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(_charmInput);
			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000.";
			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0]";
			break;
		case jetAll5:
			std::cout << "getting pT correction factor for all b-jets in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000.";
			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0]";
			break;
		case jetRecoD04:
			std::cout << "getting pT correction factor for c->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return 1.;//no correction needed
		case jetRecoSV4:
			std::cout << "getting pT correction factor for c->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(_charmInput);
			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000. && SVM[0]";
			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && SVM[0]";
			break;
		case jetRecoD05:
			std::cout << "getting pT correction factor for b->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(_beautyInput);
			//cut1+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEBPT>5000.";
			//cut2+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEDPT>5000.";
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetRecoSV5:
			std::cout << "getting pT correction factor for b->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(_beautyInput);
			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000. && SVM[0]";
			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && SVM[0]";
			break;
		default:
			return 1.;
	}

	if(!f) return 1.;
	t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return 1.;

	corr = t->GetEntries(cut1)/static_cast<double>(t->GetEntries(cut2));
	std::cout << "Correction factor: " << corr << std::endl;
	f->Close();
	return corr;
}

bool MCJets::weightMC(jetType type, TH2D* effHist) {
	if(!_inputsSet) {
		std::cout << "ERROR in MCJets::weightMC: input files not set yet" << std::endl;
		return false;
	}

	TString nameStr=typeName(type);

	TH2D *effNum(0), *effDen(0);
	if(effHist) {
		effNum = cloneTH2D(TString(effHist->GetName())+"_num",effHist);
		effDen = cloneTH2D(TString(effHist->GetName())+"_den",effHist);
		effNum->Sumw2();
		effDen->Sumw2();
	}

	TFile* fin(0);
	std::cout << "INFO : weighting " << nameStr << " MC sample" << std::endl;
	switch(type) {
		case jetAll0:
			fin = TFile::Open(_lightInput);
			break;
		case jetAll4:
		case jetRecoD04:
		case jetRecoSV4:
			fin = TFile::Open(_charmInput);
			break;
		case jetAll5:
		case jetRecoD05:
		case jetRecoSV5:
			fin = TFile::Open(_beautyInput);
			break;
		default:
			return false;
	}
	if(!fin) return 0;
	TTree* tin = dynamic_cast<TTree*>(fin->Get("T"));
	if(!tin) return 0;

	double JetPT;
	double JetTruePT;
	double JetTrueD0;
	double JetTrueDSV;
	double JetTrueBSV;

	double ZTRUEPZ, ZTRUEE, ZTRUEY;

	std::vector<double>* vSVMCor = new std::vector<double>();
	std::vector<double>* vSVN = new std::vector<double>();

	std::vector<double>* vD0M = new std::vector<double>();
	std::vector<double>* vD0PT = new std::vector<double>();
	std::vector<double>* vD0IPCHI2 = new std::vector<double>();
	std::vector<double>* vD0KWEIGHT = new std::vector<double>();
	std::vector<double>* vD0PIWEIGHT = new std::vector<double>();

	tin->SetBranchAddress("JetPT",       &JetPT);
	tin->SetBranchAddress("JetTruePT",   &JetTruePT);
	tin->SetBranchAddress("JetTRUED0",   &JetTrueD0);
	tin->SetBranchAddress("JetTRUEDSV",  &JetTrueDSV);
	tin->SetBranchAddress("JetTRUEBSV",  &JetTrueBSV);

	tin->SetBranchAddress("ZTRUEPZ",   &ZTRUEPZ);
	tin->SetBranchAddress("ZTRUEE",   &ZTRUEE);

	tin->SetBranchAddress("SVMCor",      &vSVMCor);
	tin->SetBranchAddress("SVN",         &vSVN);

	tin->SetBranchAddress("D0M",         &vD0M);
	tin->SetBranchAddress("D0PT",        &vD0PT);
	tin->SetBranchAddress("D0IPCHI2",    &vD0IPCHI2);
	tin->SetBranchAddress("D0KWEIGHT",   &vD0KWEIGHT);
	tin->SetBranchAddress("D0PIWEIGHT",  &vD0PIWEIGHT);

	//weight MC for continuous true jet pT
	if(_truePtWeights.find(type)==_truePtWeights.end()) {
		std::cout << "INFO in MCJets:weightMC: input true pT weights not set" << std::endl;
		std::cout << "                           using weight of 1" << std::endl;
		_truePtWeights[type] = new TH1D("jetTruePtWeights","",1,10e3,100e3);
		_truePtWeights[type]->SetBinContent(1,1.);
	}
	//unsigned int npt=4;
	//double* ptInputWeightBins  = new double[npt +1]{10000.,15000.,20000.,50000.,100000.};

	//TH1D* jetTruePtWeights = new TH1D("jetTruePtWeights","",npt,ptInputWeightBins);
	//if(type==jetRecoD04 || type==jetRecoSV4) {
	//	jetTruePtWeights->SetBinContent(1,1.);
	//	jetTruePtWeights->SetBinContent(2,0.6);
	//	jetTruePtWeights->SetBinContent(3,0.6);
	//	jetTruePtWeights->SetBinContent(4,0.07);
	//} else if(type==jetRecoD05 || type==jetRecoSV5) {
	//	jetTruePtWeights->SetBinContent(1,1.);
	//	jetTruePtWeights->SetBinContent(2,0.6);
	//	jetTruePtWeights->SetBinContent(3,0.5);
	//	jetTruePtWeights->SetBinContent(4,0.02);
	//} else {
	//	return false;
	//}
	//end true(pT) weights setup

	TH1D* fPtDWeights = new TH1D("fPtDWeights","",20,0.,2.0);
	TH1D* ptRecWeights = new TH1D("ptRecWeights","",90,10.e3,100.e3);

	//for D0 also weight to match pT(D0)/pT(jet) to data
	if(!_d0DataInput.IsNull() && (type==jetRecoD04 || type==jetRecoD05)) {
		TFile* fdata = TFile::Open(_d0DataInput);
		TTree* tdata = dynamic_cast<TTree*>(fdata->Get("T"));

		TH1D* fPtDData = new TH1D("fPtDData","",20,0.,2.0);
		TH1D* fPtDSim = new TH1D("fPtDSim","",20,0.,2.0);

		fPtDData->Sumw2();
		fPtDSim->Sumw2();

		tin->Draw("D0PT/JetPT>>fPtDSim","D0M>0 && JetTRUED0 && D0KPNNK>0.2 && D0PIPNNPI>0.1");
		if(type==jetRecoD04) {
			tdata->Draw("D0PT/JetPT>>fPtDData","weight4*(TMath::Abs(D0M-1865)<30.)*(D0LogIPChi2<2.5)");
		} else if(type==jetRecoD05) {
			tdata->Draw("D0PT/JetPT>>fPtDData","weight5*(TMath::Abs(D0M-1865)<30.)*(D0LogIPChi2>2.5)");
		}

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
	}
	//end f(pT) weights setup

	//for SV also weight in pT_reco(jet) to data
	if(!_svDataInput.IsNull() && (type==jetRecoSV4 || type==jetRecoSV5)) {
		TFile* fdata = TFile::Open(_svDataInput);
		TTree* tdata = dynamic_cast<TTree*>(fdata->Get("T"));

		TH1D* ptRecData = new TH1D("ptRecData","",90,10.e3,100.e3);
		TH1D* ptRecSim  = new TH1D("ptRecSim", "",90,10.e3,100.e3);

		ptRecData->Sumw2();
		ptRecSim->Sumw2();

		tdata->Draw("JetPT>>ptRecData","SVM[0]>0");
		//if(type==jetRecoSV4) {
		//	tin->Draw("JetPT>>ptRecSim","SVM[0]>0 && JetTRUEDSV");
		//} else if(type==jetRecoSV5) {
		//	tin->Draw("JetPT>>ptRecSim","SVM[0]>0 && JetTRUEBSV");
		//}
		int nEntries = tin->GetEntries();
		boost::progress_display progress( nEntries );
		for(int iEntry=0; iEntry<nEntries; ++iEntry) {
			++progress;
			tin->GetEntry(iEntry);
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

	TFile* fout = TFile::Open(fname,"RECREATE");
	TTree* tout = new TTree("T","");
	tout->SetDirectory(fout);

	double D0M(0.), D0PT(0.), D0LogIPChi2(0.), SVMCor(0.), SVN(0.), weight(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetTruePT",    &JetTruePT);
	tout->Branch("ZTRUEPZ",      &ZTRUEPZ);
	tout->Branch("ZTRUEE",       &ZTRUEE);
	tout->Branch("ZTRUEY",       &ZTRUEY);
	tout->Branch("D0M",          &D0M);
	tout->Branch("D0PT",         &D0PT);
	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
	tout->Branch("SVMCor",       &SVMCor);
	tout->Branch("SVN",          &SVN);
	tout->Branch("weight",       &weight);

	int nEntries = tin->GetEntries();
	boost::progress_display progress( nEntries );
	for(int iEntry=0; iEntry<nEntries; ++iEntry) {
		++progress;
		tin->GetEntry(iEntry);

		if(vD0M ->size()>0) {
			D0M  = vD0M->at(0);
			D0PT = vD0PT->at(0);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(0));
		} else {
			D0M  = 0.;
			D0PT = 0.;
			D0LogIPChi2 = 0.;
		}
		if(vSVMCor ->size()>0) {
			SVMCor = vSVMCor->at(0);
			SVN = vSVN->at(0);
		} else {
			SVMCor = -1.;
			SVN = 0.;
		}

		weight = _truePtWeights[type]->GetBinContent(_truePtWeights[type]->FindBin(JetTruePT));

		if(ZTRUEE>0.) ZTRUEY = 0.5*TMath::Log((ZTRUEE+ZTRUEPZ)/(ZTRUEE-ZTRUEPZ));
		else ZTRUEY = 3.;//default value for datasets without Z
		if(effDen) effDen->Fill(JetTruePT,ZTRUEY,weight);

		//only use matching jets
		switch(type) {
			case jetRecoD04:
			case jetRecoD05:
				if(!JetTrueD0) continue;
				if(D0PT < _d0ptmin || (_d0ptmax!=-1 && D0PT > _d0ptmax) ) continue;
				break;
			case jetRecoSV4:
				if(!JetTrueDSV) continue;
				if(SVMCor<0) continue;
				break;
			case jetRecoSV5:
				if(!JetTrueBSV) continue;
				if(SVMCor<0) continue;
				break;
			case jetAll0:
			case jetAll4:
			case jetAll5:
				break;
			default:
				return 0;
		}

		if(effNum) effNum->Fill(JetTruePT,ZTRUEY,weight);

		if(type==jetRecoSV4||type==jetRecoSV5) {
			if(!_svDataInput.IsNull()) weight *= ptRecWeights->GetBinContent(ptRecWeights->FindBin(JetPT));
		}
		if(type==jetRecoD04||type==jetRecoD05) {
			if(!_d0DataInput.IsNull()) weight *= fPtDWeights->GetBinContent(fPtDWeights->FindBin(D0PT/JetPT));

			//apply PID weights
			weight *= vD0KWEIGHT->at(0);
			weight *= vD0PIWEIGHT->at(0);
		}

		tout->Fill();
	}
	if(effHist) effHist->Divide(effNum,effDen);

	tout->SetDirectory(fout);
	tout->Write();
	tout->AutoSave();
	fout->Close();

	return true;
}

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
	double JetTruePT;
	double ZY(3.);
	double SVMCor;
	double D0PT;
	double weight(1.);

	t->SetBranchAddress("JetPT",       &JetPT);
	t->SetBranchAddress("JetTruePT",   &JetTruePT);
	if(_ybins) {
		t->SetBranchAddress("ZTRUEY",  &ZY);
	}

	t->SetBranchAddress("SVMCor",      &SVMCor);
	t->SetBranchAddress("D0PT",        &D0PT);
	t->SetBranchAddress("weight",      &weight);

	unsigned int nentries = t->GetEntries();
	boost::progress_display progress( nentries );
	for(unsigned int ientry=0; ientry<nentries; ++ientry) {
		++progress;
		t->GetEntry(ientry);

		//only use matching jets
		switch(type) {
			case jetRecoD04:
			case jetRecoD05:
				if(D0PT < _d0ptmin ||
				   (_d0ptmax!=-1 && D0PT > _d0ptmax) ) continue;
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

		truebin = binning->FindBin(JetTruePT);
		recobin = binning->FindBin(JetPT);

		if(recobin<=0 && truebin<=0) {
			++countSkip;
			continue;
		} else if(recobin<=0) {
			++countMiss;
			resp->Miss(JetTruePT,weight);
		} else if(recobin<=0) {
			++countFake;
			resp->Fake(JetPT,weight);
		} else {
			++countReal;
			resp->Fill(JetPT,JetTruePT,weight);
		}
	}

	f->Close();

	std::cout << countReal << "\t" << countFake << "\t" << countMiss << "\t" << countSkip << std::endl;

	_unfoldings[type][ybin] = resp;

	return resp;
}

TH1D* MCJets::unfold(TH1D* input, jetType type, uint ybin) {
	if(!input) return 0;

	TString nameStr=typeName(type);
	if(ybin!=0) {
		nameStr+="_";
		nameStr+=ybin;
	}
	std::cout << "INFO : unfolding " << nameStr << " results for file " << _name << std::endl;

	if(_d0ptmin!=5000. || _d0ptmax!=-1) {
		nameStr+="_D0";
		nameStr+=_d0ptmin;
		nameStr+="-";
		nameStr+=_d0ptmax;
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

	mcTrue->SetLineStyle(kDashed);
	mcMeas->SetMaximum(1.1*TMath::Max(mcMeas->GetMaximum(),mcTrue->GetMaximum()));

	TCanvas c;
	respMat->Draw("colz text");
	c.SaveAs(gSaveDir+"/unfoldingResponse_"+nameStr+".pdf");

	mcMeas->Draw();
	mcTrue->Draw("same");
	c.SaveAs(gSaveDir+"/unfoldingMC_"+nameStr+".pdf");

	//apply response
	RooUnfoldIds     unfolded(response, input, 2);
	TH1D* ret = dynamic_cast<TH1D*>(unfolded.Hreco());

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
	return unfolded;
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
		case jetRecoSV4:
			return "c2SV";
		case jetRecoD05:
			return "b2D0";
		case jetRecoSV5:
			return "b2SV";
		default:
			return "";
	}
}
