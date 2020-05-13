#include "SimDFitter.h"

#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "RooAbsDataStore.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimWSTool.h"
#include "RooSimultaneous.h"
#include "RooThresholdCategory.h"
#include "RooWorkspace.h"

#include "RooMultiKeysPdf.h"
#include "RooPromptShape.h"

#include <boost/progress.hpp>

#include "cloneHists.h"
#include "DatasetManager.h"
#include "outputFunctions.h"

void SimDFitter::setInputs(TString data, TString charm, TString beauty, TString eff, TString acc, bool isMC) {
	_dataFile  = data;
	_combFile  = data;
	_charmFile  = charm;
	_beautyFile  = beauty;
	_effFile   = eff;
	_accFile   = acc;
	_dataIsMC = isMC;

	_inputsSet=true;
	_effsSet = false;
	_initialised = false;
	_mcFitsDone = false;
}

void SimDFitter::setDPtRange(double ptmin, double ptmax) {
	_dptmin   = ptmin;
	_dptmax   = ptmax;
}

TString const SimDFitter::dFileName() {
	return gSaveDir+"/D0jets"+_name+".root";
}

void SimDFitter::initJetBinning(TH1D* ptBins, TH1D* yBins) {
	if(ptBins) {
		_ptBins.clear();
		for( auto bin=1; bin<=ptBins->GetNbinsX()+1; ++bin) {
			_ptBins.push_back(ptBins->GetBinLowEdge(bin));
		}
	}
	if(yBins) {
		_yBins.clear();
		for( auto bin=1; bin<=yBins->GetNbinsX()+1; ++bin) {
			_yBins.push_back(yBins->GetBinLowEdge(bin));
		}
	}
}

bool SimDFitter::init() {
	if(ws) delete ws;
	ws = new RooWorkspace("ws");
	bool status(true);

	status&=loadDatasets();
	status&=setupModel();
	_initialised = true;
	_mcFitsDone = false;

	return status;
}

bool SimDFitter::addEffs() {
	if(!_inputsSet) return false;

	auto status(true);

	//if file already exists and we're not recreating then skip this step
	if(_recreateInputs || gSystem->AccessPathName(dFileName())) {
		status &= addEffs(0);
		status &= addEffs(4);
		status &= addEffs(5);
		status &= addEffs(7);
	}

	if(status) {
		_effsSet = true;
		_initialised = false;
		_mcFitsDone = false;
	}

	return status;
}

bool SimDFitter::addEffs(int flavour) {

	std::cout << "INFO : adding efficiencies to file " << _name << std::endl;

	TFile* f2 = TFile::Open(_effFile);
	TFile* f3 = TFile::Open(_accFile);

	bool isMC = _dataIsMC;

	DatasetManager* dm = DatasetManager::getInstance();
	switch(flavour) {
		case 0:
			dm->setDataset(_combFile);
			break;
		case 4:
			dm->setDataset(_charmFile);
			isMC=true;
			break;
		case 5:
			dm->setDataset(_beautyFile);
			isMC=true;
			break;
		case 7:
			dm->setDataset(_dataFile);
			break;
		default:
			return false;
	}

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
	//the following histogram give a vertex-position-dependent correction to the two-track efficiency
	TH2D* hcor = dynamic_cast<TH2D*>(f3->Get("corr"));

	TH2D* hpid(0);
	//if(dataIsMC) hpid = dynamic_cast<TH2D*>(f2->Get("pidmceffD045"));
	//else hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	TH2D* hsel4 = dynamic_cast<TH2D*>(f2->Get("seleffD04"));
	TH2D* hsel5 = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
	//TH2D* hacc  = dynamic_cast<TH2D*>(f2->Get("acceffD045")); //superseded by RapidSim version
	//TH2D* hrec  = dynamic_cast<TH2D*>(f2->Get("receffD045"));

	if(!hacc || !hrec || !hcor || !hpid || !hsel4 || !hsel5) return false;

	//plot efficiency functions used
	TCanvas c1;
	hacc->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	hacc->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hrec->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	hrec->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hcor->GetXaxis()->SetTitle("#it{#rho}^{2}(#it{D}^{0}) [mm^{2}]");
	hcor->GetYaxis()->SetTitle("#it{z}(#it{D}^{0}) [mm]");
	hsel4->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	hsel4->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hsel5->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	hsel5->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hpid->GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	hpid->GetYaxis()->SetTitle("#it{#eta}(#it{D}^{0})");
	hacc->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0EffAcc.pdf");
	hrec->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0EffRec.pdf");
	hcor->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0EffRecCor.pdf");
	hsel4->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0EffSelPrompt.pdf");
	hsel5->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0EffSelDispl.pdf");
	hpid->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0EffPID.pdf");

	std::vector<double> *vD0M = new std::vector<double>();
	std::vector<double> *vD0IPCHI2 = new std::vector<double>();
	std::vector<double> *vD0PT = new std::vector<double>();
	std::vector<double> *vD0PX = new std::vector<double>();
	std::vector<double> *vD0PY = new std::vector<double>();
	std::vector<double> *vD0PZ = new std::vector<double>();
	std::vector<double> *vD0E = new std::vector<double>();
	std::vector<double> *vD0X = new std::vector<double>();
	std::vector<double> *vD0Y = new std::vector<double>();
	std::vector<double> *vD0Z = new std::vector<double>();
	std::vector<double> *vD0KP = new std::vector<double>();
	std::vector<double> *vD0KPT = new std::vector<double>();
	std::vector<double> *vD0KPX = new std::vector<double>();
	std::vector<double> *vD0KPY = new std::vector<double>();
	std::vector<double> *vD0KPZ = new std::vector<double>();
	std::vector<double> *vD0PIP = new std::vector<double>();
	std::vector<double> *vD0PIPT = new std::vector<double>();
	std::vector<double> *vD0PIPX = new std::vector<double>();
	std::vector<double> *vD0PIPY = new std::vector<double>();
	std::vector<double> *vD0PIPZ = new std::vector<double>();
	std::vector<double> *vD0KPNNK = new std::vector<double>();
	std::vector<double> *vD0PIPNNPI = new std::vector<double>();
	std::vector<double> *vD0KWEIGHT = new std::vector<double>();
	std::vector<double> *vD0PIWEIGHT = new std::vector<double>();
	std::vector<double> *vD0TRUEIDX = new std::vector<double>();
	std::vector<double> *vD0FROMB = new std::vector<double>();

	double JetPT;
	double JetEta;
	double JetTruePT;

	dm->setBranchAddress("D0M",           &vD0M);
	dm->setBranchAddress("D0IPCHI2",      &vD0IPCHI2);
	dm->setBranchAddress("D0PT",          &vD0PT);
	dm->setBranchAddress("D0PX",          &vD0PX);
	dm->setBranchAddress("D0PY",          &vD0PY);
	dm->setBranchAddress("D0PZ",          &vD0PZ);
	dm->setBranchAddress("D0E",           &vD0E);
	dm->setBranchAddress("D0X",           &vD0X);
	dm->setBranchAddress("D0Y",           &vD0Y);
	dm->setBranchAddress("D0Z",           &vD0Z);
	dm->setBranchAddress("D0KP",          &vD0KP);
	dm->setBranchAddress("D0KPT",         &vD0KPT);
	dm->setBranchAddress("D0KPX",         &vD0KPX);
	dm->setBranchAddress("D0KPY",         &vD0KPY);
	dm->setBranchAddress("D0KPZ",         &vD0KPZ);
	dm->setBranchAddress("D0PIP",         &vD0PIP);
	dm->setBranchAddress("D0PIPT",        &vD0PIPT);
	dm->setBranchAddress("D0PIPX",        &vD0PIPX);
	dm->setBranchAddress("D0PIPY",        &vD0PIPY);
	dm->setBranchAddress("D0PIPZ",        &vD0PIPZ);
	dm->setBranchAddress("D0KPNNK",       &vD0KPNNK);
	dm->setBranchAddress("D0PIPNNPI",     &vD0PIPNNPI);
	dm->setBranchAddress("D0KWEIGHT",     &vD0KWEIGHT);
	dm->setBranchAddress("D0PIWEIGHT",    &vD0PIWEIGHT);

	dm->setBranchAddress("JetPT",         &JetPT);
	dm->setBranchAddress("JetEta",        &JetEta);
	if(isMC) {
		dm->setBranchAddress("JetTruePT",     &JetTruePT);
		dm->setBranchAddress("D0TRUEIDX",     &vD0TRUEIDX);
		dm->setBranchAddress("D0FROMB",       &vD0FROMB);
	}

	unsigned int nentries0 = dm->getEntries();

	TFile* fout = TFile::Open(dFileName(),"UPDATE");
	TString tname = "T";
	tname+=flavour;
	TTree* tout = new TTree(tname,"");

	double D0M(0.), D0PT(0.), D0Eta(0.), D0LogIPChi2(0.), weight4(0.), weight5(0.);
	//TODO//double D0RhoSq(0.), D0Z(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetEta",       &JetEta);
	tout->Branch("D0M",          &D0M);
	tout->Branch("D0PT",         &D0PT);
	tout->Branch("D0Eta",        &D0Eta);
	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
	tout->Branch("weight4",      &weight4);
	tout->Branch("weight5",      &weight5);

	tout->Branch("JetTruePT", &JetTruePT);

	//first make non-vector tree for fits
	boost::progress_display progress( nentries0 );
	//for(unsigned int ientry=0; ientry<nentries0; ++ientry)
	while(dm->getNext()) {
		++progress;
		//t0->GetEntry(ientry);

		for(unsigned int s=0; s<vD0M->size(); ++s) {
			//first check in our tight acceptance
			TVector3 D0P (vD0PX->at(s)  ,vD0PY->at(s)  ,vD0PZ->at(s));
			TVector3 D0P0(vD0KPX->at(s) ,vD0KPY->at(s) ,vD0KPZ->at(s));
			TVector3 D0P1(vD0PIPX->at(s),vD0PIPY->at(s),vD0PIPZ->at(s));

			if(flavour==7 && isMC && _truthMatchData!=truthMatchType::truthMatchOff) {
				if(_truthMatchData==truthMatchType::truthMatchPrompt) {
					if(vD0TRUEIDX->at(s)==-1 || vD0FROMB->at(s)==1) continue;
				}
				if(_truthMatchData==truthMatchType::truthMatchDispl) {
					if(vD0TRUEIDX->at(s)==-1 || vD0FROMB->at(s)==0) continue;
				}
			}

			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;
			if((flavour==4 || flavour==5) && vD0TRUEIDX->at(s)==-1) continue;
			if(flavour==5 && vD0FROMB->at(s)==0) continue;
			if(flavour==4 && vD0FROMB->at(s)==1) continue;
			if(flavour==0 && vD0M->at(s)>1844. && vD0M->at(s)<1884.) continue;

			//check PID cuts
//TODO//20200420			if(!isMC && vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			if(vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			//if(!dataIsMC && vD0PIPNNPI->at(s)<0.1 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //pion PID turned off

			D0M         = vD0M->at(s);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
			D0PT        = D0P.Pt();
			D0Eta       = D0P.Eta();
			//TODO//D0RhoSq     = vD0X->at(s)*vD0X->at(s) + vD0Y->at(s)*vD0Y->at(s);
			//TODO//D0Z         = vD0Z->at(s);

			//get efficiency corrections
			double effacc(0.), effrec(0.), effcor(1.), effsel4(0.), effsel5(0.), effpid(0.);

			if(D0PT>=100000.) D0PT=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(D0PT/1000.,D0Eta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(D0PT/1000.,D0Eta));
			//TODO//if(_useRhoZCorr) effcor = hcor ->GetBinContent(hcor ->FindBin(D0RhoSq   ,D0Z  ));
			effsel4=      hsel4->GetBinContent(hsel4->FindBin(D0PT      ,D0Eta));
			effsel5=      hsel5->GetBinContent(hsel5->FindBin(D0PT      ,D0Eta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(D0PT      ,D0Eta));
//TODO//20200420			if(isMC || D0P0.Pt()>=25000. || D0P0.Mag()>=500000.) effpid=1.0; //no PID cut used in these cases
			if(D0P0.Pt()>=25000. || D0P0.Mag()>=500000.) effpid=1.0; //no PID cut used in these cases

			double eff4=effacc*effrec*effcor*effsel4*effpid;
			double eff5=effacc*effrec*effcor*effsel5*effpid;

			if(eff4<0.01 || eff5<0.01) {
//TODO//				std::cout << D0PT << "\t" << D0Eta << "\t" << effacc << "\t" << effrec << "\t" << effcor << "\t" << effsel4 << "\t" << effsel5 << "\t" << effpid << std::endl;
//TODO//				if(effcor==0) std::cout << D0RhoSq << "\t" << D0Z << std::endl;//TODO
				continue;
			}

			weight4 = 1./eff4;
			weight5 = 1./eff5;

			if(isMC) {
				//put the PID efficiency into the weights for MC
				//TODO//testWithPIDWeightsOffForTrueD0Sample//if(D0P0.Pt()<25000. && D0P0.Mag()<500000.) {
				//TODO//testWithPIDWeightsOffForTrueD0Sample//	weight4 *= vD0KWEIGHT->at(s);
				//TODO//testWithPIDWeightsOffForTrueD0Sample//	weight5 *= vD0KWEIGHT->at(s);
				//TODO//testWithPIDWeightsOffForTrueD0Sample//}
				//currently no pion PID
				//if(D0P1.Pt()<25000. && D0P1.Mag()<500000.) {
				//	weight4 *= vD0PIWEIGHT->at(s);
				//	weight5 *= vD0PIWEIGHT->at(s);
				//}
				//scale weights for roughly continuous jet true pT
				//TODO//if(JetTruePT>50000.) {//TODO added weights back in
				//TODO//	weight4*=0.03;
				//TODO//	weight5*=0.02;
				//TODO//} else if(JetTruePT>20000.) {
				//TODO//	weight4*=0.6;
				//TODO//	weight5*=0.4;
				//TODO//} else if(JetTruePT>15000.) {
				//TODO//	weight4*=0.6;
				//TODO//	weight5*=0.5;
				//TODO//}
				//scale weights for roughly continuous jet true pT
				//TODO//breaks eff if comparing against unweighted MC//if(JetTruePT>50000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.007;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.007;
				//TODO//breaks eff if comparing against unweighted MC//} else if(JetTruePT>20000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.10;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.10;
				//TODO//breaks eff if comparing against unweighted MC//} else if(JetTruePT>15000.) {
				//TODO//breaks eff if comparing against unweighted MC//	weight4*=0.25;
				//TODO//breaks eff if comparing against unweighted MC//	weight5*=0.25;
				//TODO//breaks eff if comparing against unweighted MC//}
			}       //TODO//breaks eff if comparing against unweighted MC//
                                //TODO//breaks eff if comparing against unweighted MC//
			tout->Fill();
			//TODO//20200420-keepAllD0s//break;//only keep one D0 candidate per entry
		}
	}

	dm->reset();

	tout->AutoSave();
	fout->Close();

	return true;
}

bool SimDFitter::testEffs(int flavour) {
	if(!_dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << _name << std::endl;
	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset(_dataFile);

	TFile* f2 = TFile::Open(_effFile);
	TFile* f3 = TFile::Open(_accFile);

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
	TH2D* hcor = dynamic_cast<TH2D*>(f3->Get("corr"));

	TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	TH2D* hsel(0);
	if(_flavour==5) hsel = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
	else hsel = dynamic_cast<TH2D*>(f2->Get("seleffD04"));
	//TH2D* hacc  = dynamic_cast<TH2D*>(f2->Get("acceffD045"));
	//TH2D* hrec  = dynamic_cast<TH2D*>(f2->Get("receffD045"));

	if(!hacc || !hrec || !hcor || !hpid || !hsel) return false;

	//add over/underflow to the pT binning
	std::vector<double> ptBinsWithOF = {0.};
	for(uint i=0; i<_ptBins.size(); ++i) ptBinsWithOF.push_back(_ptBins.at(i));
	ptBinsWithOF.push_back(110e3);
	for(uint i=0; i<ptBinsWithOF.size(); ++i) std::cout << ptBinsWithOF.at(i) << std::endl;//TODO

	//yield histograms - 0=all; 1=passAcc; 2=passRec; 3=passSel; 4=passPID.
	//tru=truth D0; fit=eff-corrected reco D0; fnd=eff-corrected truth-matched reco D0.
	std::vector<TString> stages;
	stages.push_back("all");
	stages.push_back("geometric");
	stages.push_back("reconstruction");
	stages.push_back("selection");
	stages.push_back("PID");
	std::vector<TH1D*> truD0truePT;//the true distribution at this stage
	truD0truePT.push_back(cloneTH1D("tru0truept", ptBinsWithOF));
	truD0truePT.push_back(cloneTH1D("tru1truept", ptBinsWithOF));
	truD0truePT.push_back(cloneTH1D("tru2truept", ptBinsWithOF));
	truD0truePT.push_back(cloneTH1D("tru3truept", ptBinsWithOF));
	truD0truePT.push_back(cloneTH1D("tru4truept", ptBinsWithOF));
	std::vector<TH1D*> effD0truePT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effD0truePT.push_back(cloneTH1D("eff0truept", ptBinsWithOF));
	effD0truePT.push_back(cloneTH1D("eff1truept", ptBinsWithOF));
	effD0truePT.push_back(cloneTH1D("eff2truept", ptBinsWithOF));
	effD0truePT.push_back(cloneTH1D("eff3truept", ptBinsWithOF));
	effD0truePT.push_back(cloneTH1D("eff4truept", ptBinsWithOF));
	std::vector<TH1D*> fitD0truePT;//sideband-subtracted distribution at this stage
	fitD0truePT.push_back(cloneTH1D("fit0truept", ptBinsWithOF));
	fitD0truePT.push_back(cloneTH1D("fit1truept", ptBinsWithOF));
	fitD0truePT.push_back(cloneTH1D("fit2truept", ptBinsWithOF));
	fitD0truePT.push_back(cloneTH1D("fit3truept", ptBinsWithOF));
	fitD0truePT.push_back(cloneTH1D("fit4truept", ptBinsWithOF));
	std::vector<TH1D*> fndD0truePT;//truth-matched distribution at this stage
	fndD0truePT.push_back(cloneTH1D("fnd0truept", ptBinsWithOF));
	fndD0truePT.push_back(cloneTH1D("fnd1truept", ptBinsWithOF));
	fndD0truePT.push_back(cloneTH1D("fnd2truept", ptBinsWithOF));
	fndD0truePT.push_back(cloneTH1D("fnd3truept", ptBinsWithOF));
	fndD0truePT.push_back(cloneTH1D("fnd4truept", ptBinsWithOF));
	//binned in true jet pT
	std::vector<TH1D*> truD0recoPT;
	truD0recoPT.push_back(cloneTH1D("tru0recopt", ptBinsWithOF));
	truD0recoPT.push_back(cloneTH1D("tru1recopt", ptBinsWithOF));
	truD0recoPT.push_back(cloneTH1D("tru2recopt", ptBinsWithOF));
	truD0recoPT.push_back(cloneTH1D("tru3recopt", ptBinsWithOF));
	truD0recoPT.push_back(cloneTH1D("tru4recopt", ptBinsWithOF));
	std::vector<TH1D*> effD0recoPT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effD0recoPT.push_back(cloneTH1D("eff0recopt", ptBinsWithOF));
	effD0recoPT.push_back(cloneTH1D("eff1recopt", ptBinsWithOF));
	effD0recoPT.push_back(cloneTH1D("eff2recopt", ptBinsWithOF));
	effD0recoPT.push_back(cloneTH1D("eff3recopt", ptBinsWithOF));
	effD0recoPT.push_back(cloneTH1D("eff4recopt", ptBinsWithOF));
	std::vector<TH1D*> fitD0recoPT;
	fitD0recoPT.push_back(cloneTH1D("fit0recopt", ptBinsWithOF));
	fitD0recoPT.push_back(cloneTH1D("fit1recopt", ptBinsWithOF));
	fitD0recoPT.push_back(cloneTH1D("fit2recopt", ptBinsWithOF));
	fitD0recoPT.push_back(cloneTH1D("fit3recopt", ptBinsWithOF));
	fitD0recoPT.push_back(cloneTH1D("fit4recopt", ptBinsWithOF));
	std::vector<TH1D*> fndD0recoPT;
	fndD0recoPT.push_back(cloneTH1D("fnd0recopt", ptBinsWithOF));
	fndD0recoPT.push_back(cloneTH1D("fnd1recopt", ptBinsWithOF));
	fndD0recoPT.push_back(cloneTH1D("fnd2recopt", ptBinsWithOF));
	fndD0recoPT.push_back(cloneTH1D("fnd3recopt", ptBinsWithOF));
	fndD0recoPT.push_back(cloneTH1D("fnd4recopt", ptBinsWithOF));
	//binned in both reco and true jet pT
	std::vector<TH2D*> truD02DPT;
	truD02DPT.push_back(cloneTH2D("tru02Dpt", ptBinsWithOF, ptBinsWithOF));
	truD02DPT.push_back(cloneTH2D("tru12Dpt", ptBinsWithOF, ptBinsWithOF));
	truD02DPT.push_back(cloneTH2D("tru22Dpt", ptBinsWithOF, ptBinsWithOF));
	truD02DPT.push_back(cloneTH2D("tru32Dpt", ptBinsWithOF, ptBinsWithOF));
	truD02DPT.push_back(cloneTH2D("tru42Dpt", ptBinsWithOF, ptBinsWithOF));
	std::vector<TH2D*> effD02DPT;
	effD02DPT.push_back(cloneTH2D("eff02Dpt", ptBinsWithOF, ptBinsWithOF));
	effD02DPT.push_back(cloneTH2D("eff12Dpt", ptBinsWithOF, ptBinsWithOF));
	effD02DPT.push_back(cloneTH2D("eff22Dpt", ptBinsWithOF, ptBinsWithOF));
	effD02DPT.push_back(cloneTH2D("eff32Dpt", ptBinsWithOF, ptBinsWithOF));
	effD02DPT.push_back(cloneTH2D("eff42Dpt", ptBinsWithOF, ptBinsWithOF));
	std::vector<TH2D*> fitD02DPT;
	fitD02DPT.push_back(cloneTH2D("fit02Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD02DPT.push_back(cloneTH2D("fit12Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD02DPT.push_back(cloneTH2D("fit22Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD02DPT.push_back(cloneTH2D("fit32Dpt", ptBinsWithOF, ptBinsWithOF));
	fitD02DPT.push_back(cloneTH2D("fit42Dpt", ptBinsWithOF, ptBinsWithOF));
	std::vector<TH2D*> fndD02DPT;
	fndD02DPT.push_back(cloneTH2D("fnd02Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD02DPT.push_back(cloneTH2D("fnd12Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD02DPT.push_back(cloneTH2D("fnd22Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD02DPT.push_back(cloneTH2D("fnd32Dpt", ptBinsWithOF, ptBinsWithOF));
	fndD02DPT.push_back(cloneTH2D("fnd42Dpt", ptBinsWithOF, ptBinsWithOF));

	//entries in this vector are true and reco for each stage
	std::vector<TH2D*> d0Kinematics;
	d0Kinematics.push_back(cloneTH2D("d0KineTrue0", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineReco0", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue1", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineReco1", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue2", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineReco2", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue3", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineReco3", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue4", hsel));
	d0Kinematics.push_back(cloneTH2D("d0KineReco4", hsel));
	for(unsigned int i=0; i<d0Kinematics.size(); ++i) {
		d0Kinematics[i]->Sumw2();
	}
	std::vector<TH1D*> d0KinematicsPulls;
	d0KinematicsPulls.push_back(new TH1D("pulls1", "", 30, -1., 1.));
	d0KinematicsPulls.push_back(new TH1D("pulls2", "", 30, -1., 1.));
	d0KinematicsPulls.push_back(new TH1D("pulls3", "", 30, -1., 1.));
	d0KinematicsPulls.push_back(new TH1D("pulls4", "", 30, -1., 1.));

	//also divide distribution into bins of true and reco pT
	std::vector<TH2D*> d0KinematicsTRBinned;
	std::vector<TH2D*> d0KinematicsTRBinned_recD_tru;
	std::vector<TH2D*> d0KinematicsTRBinned_recD_rec;
	for(uint i=0; i<ptBinsWithOF.size(); ++i) {
		for(uint j=0; j<ptBinsWithOF.size()-1; ++j) {
			TString hName = TString::Format("d0KineJPTBin%d%d",i,j);
			TString hTitle = TString::Format("true%d reco%d",j,i);
			if(j==5) {
				hTitle.Form("reco%d",i);
			}
			if(i==6) {
				hTitle.Form("true%d",j);
			}
			if(j==5 && i==6) {
				hTitle = "all";
			}
			d0KinematicsTRBinned.push_back(cloneTH2D(hName, hsel));
			d0KinematicsTRBinned.back()->Sumw2();
			d0KinematicsTRBinned.back()->SetTitle(hTitle);
			d0KinematicsTRBinned.back()->GetXaxis()->SetTitle("pT(D)");
			d0KinematicsTRBinned.back()->GetYaxis()->SetTitle("#eta(D)");
			
			d0KinematicsTRBinned_recD_tru.push_back(cloneTH2D(hName+"_recD_tru", hsel));
			d0KinematicsTRBinned_recD_tru.back()->Sumw2();
			d0KinematicsTRBinned_recD_tru.back()->SetTitle(hTitle);
			d0KinematicsTRBinned_recD_tru.back()->GetXaxis()->SetTitle("pT(D)");
			d0KinematicsTRBinned_recD_tru.back()->GetYaxis()->SetTitle("#eta(D)");
			
			d0KinematicsTRBinned_recD_rec.push_back(cloneTH2D(hName+"_recD_rec", hsel));
			d0KinematicsTRBinned_recD_rec.back()->Sumw2();
			d0KinematicsTRBinned_recD_rec.back()->SetTitle(hTitle);
			d0KinematicsTRBinned_recD_rec.back()->GetXaxis()->SetTitle("pT(D)");
			d0KinematicsTRBinned_recD_rec.back()->GetYaxis()->SetTitle("#eta(D)");
		}
	}

	for(int i=0; i<5; ++i) {
		truD0truePT.at(i)->Sumw2();
		effD0truePT.at(i)->Sumw2();
		fitD0truePT.at(i)->Sumw2();
		fndD0truePT.at(i)->Sumw2();
		truD0recoPT.at(i)->Sumw2();
		effD0recoPT.at(i)->Sumw2();
		fitD0recoPT.at(i)->Sumw2();
		fndD0recoPT.at(i)->Sumw2();
		truD02DPT.at(i)->Sumw2();
		effD02DPT.at(i)->Sumw2();
		fitD02DPT.at(i)->Sumw2();
		fndD02DPT.at(i)->Sumw2();
		truD02DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		effD02DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		fitD02DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		fndD02DPT.at(i)->GetXaxis()->SetTitle("true jet pT");
		truD02DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
		effD02DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
		fitD02DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
		fndD02DPT.at(i)->GetYaxis()->SetTitle("reco jet pT");
	}

	std::vector<double> *D0M = new std::vector<double>();
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
	std::vector<double> *D0KWEIGHT = new std::vector<double>();
	std::vector<double> *D0PIWEIGHT = new std::vector<double>();
	std::vector<double> *D0TRUEIDX = new std::vector<double>();
	std::vector<double> *D0TRUETRK0 = new std::vector<double>();
	std::vector<double> *D0TRUETRK1 = new std::vector<double>();

	std::vector<double> *TRUEDID = new std::vector<double>();
	std::vector<double> *TRUEDPX = new std::vector<double>();
	std::vector<double> *TRUEDPY = new std::vector<double>();
	std::vector<double> *TRUEDPZ = new std::vector<double>();
	std::vector<double> *TRUEDE = new std::vector<double>();
	std::vector<double> *TRUEDX = new std::vector<double>();
	std::vector<double> *TRUEDY = new std::vector<double>();
	std::vector<double> *TRUEDZ = new std::vector<double>();
	std::vector<double> *TRUEDFROMB = new std::vector<double>();
	std::vector<double> *TRUEDTRK0P = new std::vector<double>();
	std::vector<double> *TRUEDTRK0PT = new std::vector<double>();
	std::vector<double> *TRUEDTRK0INACC = new std::vector<double>();
	std::vector<double> *TRUEDTRK0RECO = new std::vector<double>();
	std::vector<double> *TRUEDTRK1P = new std::vector<double>();
	std::vector<double> *TRUEDTRK1PT = new std::vector<double>();
	std::vector<double> *TRUEDTRK1INACC = new std::vector<double>();
	std::vector<double> *TRUEDTRK1RECO = new std::vector<double>();
	std::vector<double> *TRUEDTRK0IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK1IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK2IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK3IDX = new std::vector<double>();
	std::vector<double> *TRUEDTRK0ID = new std::vector<double>();
	std::vector<double> *TRUEDTRK1ID = new std::vector<double>();
	std::vector<double> *TRUEDTRK2ID = new std::vector<double>();
	std::vector<double> *TRUEDTRK3ID = new std::vector<double>();

	double JetPT;
	double JetEta;
	double JetTruePT;

	dm->setBranchAddress("D0M",           &D0M);
	dm->setBranchAddress("D0IPCHI2",      &D0IPCHI2);
	dm->setBranchAddress("D0PT",          &D0PT);
	dm->setBranchAddress("D0PX",          &D0PX);
	dm->setBranchAddress("D0PY",          &D0PY);
	dm->setBranchAddress("D0PZ",          &D0PZ);
	dm->setBranchAddress("D0E",           &D0E);
	dm->setBranchAddress("D0X",           &D0X);
	dm->setBranchAddress("D0Y",           &D0Y);
	dm->setBranchAddress("D0Z",           &D0Z);
	dm->setBranchAddress("D0KP",          &D0KP);
	dm->setBranchAddress("D0KPT",         &D0KPT);
	dm->setBranchAddress("D0KPX",         &D0KPX);
	dm->setBranchAddress("D0KPY",         &D0KPY);
	dm->setBranchAddress("D0KPZ",         &D0KPZ);
	dm->setBranchAddress("D0PIP",         &D0PIP);
	dm->setBranchAddress("D0PIPT",        &D0PIPT);
	dm->setBranchAddress("D0PIPX",        &D0PIPX);
	dm->setBranchAddress("D0PIPY",        &D0PIPY);
	dm->setBranchAddress("D0PIPZ",        &D0PIPZ);
	dm->setBranchAddress("D0KWEIGHT",     &D0KWEIGHT);
	dm->setBranchAddress("D0PIWEIGHT",    &D0PIWEIGHT);
	dm->setBranchAddress("D0TRUEIDX",     &D0TRUEIDX);
	dm->setBranchAddress("D0TRUETRK0",    &D0TRUETRK0);
	dm->setBranchAddress("D0TRUETRK1",    &D0TRUETRK1);

	dm->setBranchAddress("TRUEDID",       &TRUEDID);
	dm->setBranchAddress("TRUEDPX",       &TRUEDPX);
	dm->setBranchAddress("TRUEDPY",       &TRUEDPY);
	dm->setBranchAddress("TRUEDPZ",       &TRUEDPZ);
	dm->setBranchAddress("TRUEDE",        &TRUEDE);
	dm->setBranchAddress("TRUEDX",        &TRUEDX);
	dm->setBranchAddress("TRUEDY",        &TRUEDY);
	dm->setBranchAddress("TRUEDZ",        &TRUEDZ);
	dm->setBranchAddress("TRUEDFROMB",    &TRUEDFROMB);
	dm->setBranchAddress("TRUEDTRK0P",    &TRUEDTRK0P);
	dm->setBranchAddress("TRUEDTRK0PT",   &TRUEDTRK0PT);
	dm->setBranchAddress("TRUEDTRK0INACC",&TRUEDTRK0INACC);
	dm->setBranchAddress("TRUEDTRK0RECO", &TRUEDTRK0RECO);
	dm->setBranchAddress("TRUEDTRK1P",    &TRUEDTRK1P);
	dm->setBranchAddress("TRUEDTRK1PT",   &TRUEDTRK1PT);
	dm->setBranchAddress("TRUEDTRK1INACC",&TRUEDTRK1INACC);
	dm->setBranchAddress("TRUEDTRK1RECO", &TRUEDTRK1RECO);
	dm->setBranchAddress("TRUEDTRK0IDX",  &TRUEDTRK0IDX);
	dm->setBranchAddress("TRUEDTRK1IDX",  &TRUEDTRK1IDX);
	dm->setBranchAddress("TRUEDTRK2IDX",  &TRUEDTRK2IDX);
	dm->setBranchAddress("TRUEDTRK3IDX",  &TRUEDTRK3IDX);
	dm->setBranchAddress("TRUEDTRK0ID",  &TRUEDTRK0ID);
	dm->setBranchAddress("TRUEDTRK1ID",  &TRUEDTRK1ID);
	dm->setBranchAddress("TRUEDTRK2ID",  &TRUEDTRK2ID);
	dm->setBranchAddress("TRUEDTRK3ID",  &TRUEDTRK3ID);

	dm->setBranchAddress("JetPT",         &JetPT);
	dm->setBranchAddress("JetEta",        &JetEta);
	dm->setBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = dm->getEntries();

	double dpt(0.), deta(0.);

	boost::progress_display progress( nentries0 );
	//for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
	while(dm->getNext()) {
		++progress;
		//t0->GetEntry(ientry);

		std::vector<int> trueD0s, foundD0s;

		//put anything out of range into the overflow bin
		if(JetPT>_ptBins.back()) JetPT = _ptBins.back();
		if(JetTruePT>_ptBins.back()) JetTruePT = _ptBins.back();

		for(unsigned int d=0; d<TRUEDID->size(); ++d) {
			//only use D0->Kpi with Pt>5GeV
			if(TMath::Abs(TRUEDID->at(d))!=421) continue;
			if(TRUEDTRK0IDX->at(d)==-1) continue;
			if(TRUEDTRK2IDX->at(d)!=-1) continue;
			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;

			TVector3 TRUEDP (TRUEDPX->at(d)  ,TRUEDPY->at(d)  ,TRUEDPZ->at(d));

			dpt        = TRUEDP.Pt();
			deta       = TRUEDP.Eta();
			double rhoSq = TRUEDX->at(d)*TRUEDX->at(d) + TRUEDY->at(d)*TRUEDY->at(d);
			double z = TRUEDZ->at(d);
			if(rhoSq>=100) rhoSq=99.9;//TODO

			//get efficiency corrections
			double effacc(0.), effrec(0.), effsel(0.), effpid(0.);

			if(dpt>=100000.) dpt=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt/1000.,deta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt/1000.,deta));
			//apply correction to the reco efficiency
			if(false) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));

			truD0recoPT.at(0)->Fill(JetPT);
			truD0truePT.at(0)->Fill(JetTruePT);
			truD02DPT.at(0)->Fill(JetTruePT,JetPT);
			d0Kinematics[0]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			//fill kinematics separately for each true/reco jet pT bin
			int ibin = truD0recoPT.at(0)->FindBin(JetTruePT) - 1;
			int jbin = truD0recoPT.at(0)->FindBin(JetPT) - 1;
			d0KinematicsTRBinned.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());
			d0KinematicsTRBinned.at(5 + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//row sums
			d0KinematicsTRBinned.at(ibin + 36)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//col sums
			d0KinematicsTRBinned.at(41)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//all

			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;

			truD0recoPT.at(1)->Fill(JetPT);
			truD0truePT.at(1)->Fill(JetTruePT);
			truD02DPT.at(1)->Fill(JetTruePT,JetPT);
			if(effacc>0.) {
				effD0recoPT.at(0)->Fill(JetPT,1./effacc);
				effD0truePT.at(0)->Fill(JetTruePT,1./effacc);
				effD02DPT.at(0)->Fill(JetTruePT,JetPT,1./effacc);
			} else {
				std::cout << "ACC=" << effacc << " at " << dpt << "," << deta << std::endl;
			}
			d0Kinematics[2]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			truD0recoPT.at(2)->Fill(JetPT);
			truD0truePT.at(2)->Fill(JetTruePT);
			truD02DPT.at(2)->Fill(JetTruePT,JetPT);
			if(effrec>0.) {
				effD0recoPT.at(1)->Fill(JetPT,1./effrec);
				effD0truePT.at(1)->Fill(JetTruePT,1./effrec);
				effD02DPT.at(1)->Fill(JetTruePT,JetPT,1./effrec);
			} else {
				std::cout << "REC=" << effrec << " at " << dpt << "," << deta << std::endl;
			}
			d0Kinematics[4]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			for(unsigned int s=0; s<D0TRUEIDX->size(); ++s) {
				if(D0TRUEIDX->at(s)==d) {
					TVector3 D0P (D0PX->at(s)  ,D0PY->at(s)  ,D0PZ->at(s));
					TVector3 D0P0(D0KPX->at(s) ,D0KPY->at(s) ,D0KPZ->at(s));
					TVector3 D0P1(D0PIPX->at(s),D0PIPY->at(s),D0PIPZ->at(s));

					if(!(D0P0.Eta()>2.0 && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
					if(!(D0P1.Eta()>2.0 && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
					if(!(D0P.Pt()>_dptmin)) continue;
					if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;
					//if(!(D0P.Eta()>2.5&&D0P.Eta()<4.0)) continue;//TODO test restricting the eta(D0) range (turn on in two places 1/2)

					truD0recoPT.at(3)->Fill(JetPT);
					truD0truePT.at(3)->Fill(JetTruePT);
					truD02DPT.at(3)->Fill(JetTruePT,JetPT);
					if(effsel>0.) {
						effD0recoPT.at(2)->Fill(JetPT,1./effsel);
						effD0truePT.at(2)->Fill(JetTruePT,1./effsel);
						effD02DPT.at(2)->Fill(JetTruePT,JetPT,1./effsel);
					} else {
						std::cout << "SEL=" << effsel << " at " << dpt << "," << deta << std::endl;
					}
					d0Kinematics[6]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

					//use weights for PID
					double weight = D0KWEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
					if(D0P0.Pt()>25000. || D0P0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT
					truD0recoPT.at(4)->Fill(JetPT,weight);
					truD0truePT.at(4)->Fill(JetTruePT,weight);
					truD02DPT.at(4)->Fill(JetTruePT,JetPT);
					if(effpid>0.) {
						effD0recoPT.at(3)->Fill(JetPT,weight/effpid);
						effD0truePT.at(3)->Fill(JetTruePT,weight/effpid);
						effD02DPT.at(3)->Fill(JetTruePT,JetPT,1./effpid);
					} else {
						std::cout << "PID=" << effpid << " at " << dpt << "," << deta << std::endl;
					}
					d0Kinematics[8]->Fill(TRUEDP.Pt(),TRUEDP.Eta(),weight);

					trueD0s.push_back(d*100+s);

					//fill kinematics separately for each true/reco jet pT bin
					d0KinematicsTRBinned_recD_tru.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());
					d0KinematicsTRBinned_recD_tru.at(5 + jbin*(ptBinsWithOF.size()-1))->Fill(TRUEDP.Pt(),TRUEDP.Eta());//row sums
					d0KinematicsTRBinned_recD_tru.at(ibin + 36)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//col sums
					d0KinematicsTRBinned_recD_tru.at(41)->Fill(TRUEDP.Pt(),TRUEDP.Eta());//all

					//fill kinematics separately for each true/reco jet pT bin
					d0KinematicsTRBinned_recD_rec.at(ibin + jbin*(ptBinsWithOF.size()-1))->Fill(D0P.Pt(),D0P.Eta());
					d0KinematicsTRBinned_recD_rec.at(5 + jbin*(ptBinsWithOF.size()-1))->Fill(D0P.Pt(),D0P.Eta());//row sums
					d0KinematicsTRBinned_recD_rec.at(ibin + 36)->Fill(D0P.Pt(),D0P.Eta());//col sums
					d0KinematicsTRBinned_recD_rec.at(41)->Fill(D0P.Pt(),D0P.Eta());//all

					break;
				}
			}
		}
		for(unsigned int s=0; s<D0M->size(); ++s) {
			//first check in our tight acceptance
			TVector3 D0P (D0PX->at(s)  ,D0PY->at(s)  ,D0PZ->at(s));
			TVector3 D0P0(D0KPX->at(s) ,D0KPY->at(s) ,D0KPZ->at(s));
			TVector3 D0P1(D0PIPX->at(s),D0PIPY->at(s),D0PIPZ->at(s));
			double rhoSq = D0X->at(s)*D0X->at(s) + D0Y->at(s)*D0Y->at(s);
			double z = D0Z->at(s);
			if(rhoSq>=100) rhoSq=99.9;//TODO

			if(!(D0P0.Eta()>2.0 && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2.0 && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;
			//if(!(D0P.Eta()>2.2&&D0P.Eta()<4.0)) continue;//TODO test restrcting the eta(D0) range (turn on in two places 2/2)

			//check PID cuts
			double weight = D0KWEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
			if(D0P0.Pt()>25000. || D0P0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT

			dpt        = D0P.Pt();
			deta       = D0P.Eta();

			//get efficiency corrections
			double effacc(0.), effrec(0.), effsel(0.), effpid(0.);

			if(dpt>=100000.) dpt=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(dpt/1000.,deta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(dpt/1000.,deta));
			//apply correction to the reco efficiency
			if(false) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));

			double eff=effacc*effrec*effsel*effpid;

			if(eff<0.01) {
				//TODO//std::cout << dpt << "\t" << deta << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
				//if(D0TRUEIDX->at(s)>-1) {
				//	int d = D0TRUEIDX->at(s);
				//	if(TRUEDTRK0IDX->at(d)!=-1 && TRUEDTRK2IDX->at(d)==-1 && TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)>=5000.*5000.) {
				//		if(TRUEDTRK0INACC->at(d)==1 && TRUEDTRK1INACC->at(d)==1) {
				//			if(TRUEDTRK0RECO->at(d)==1 && TRUEDTRK1RECO->at(d)==1) {
				//				std::cout << "LOST TRUE D0" << std::endl;
				//std::cout << dpt << "\t" << deta << "\t" << eff << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
				//std::cout << rhoSq << "\t" << z << "\t" << hrec ->GetBinContent(hrec ->FindBin(dpt/1000.,deta)) << "\t" << hcor ->GetBinContent(hcor ->FindBin(rhoSq,z)) << std::endl;
				//			}
				//		}
				//	}
				//}
				continue;
			}

			double sbSubSwitch(0.);
			if(TMath::Abs(D0M->at(s)-1865)<40.) sbSubSwitch = 1.;
			else if(TMath::Abs(D0M->at(s)-1865)<80.) sbSubSwitch = -1.;

			fitD0recoPT.at(4)->Fill(JetPT,weight*sbSubSwitch);
			fitD0recoPT.at(3)->Fill(JetPT,weight*sbSubSwitch/effpid);
			fitD0recoPT.at(2)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel));
			fitD0recoPT.at(1)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel*effrec));
			fitD0recoPT.at(0)->Fill(JetPT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			fitD0truePT.at(4)->Fill(JetTruePT,weight*sbSubSwitch);
			fitD0truePT.at(3)->Fill(JetTruePT,weight*sbSubSwitch/effpid);
			fitD0truePT.at(2)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel));
			fitD0truePT.at(1)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel*effrec));
			fitD0truePT.at(0)->Fill(JetTruePT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			fitD02DPT.at(4)->Fill(JetTruePT,JetPT,weight*sbSubSwitch);
			fitD02DPT.at(3)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/effpid);
			fitD02DPT.at(2)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel));
			fitD02DPT.at(1)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel*effrec));
			fitD02DPT.at(0)->Fill(JetTruePT,JetPT,weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			d0Kinematics[9]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch);
			d0Kinematics[7]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch/(effpid));
			d0Kinematics[5]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch/(effpid*effsel));
			d0Kinematics[3]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch/(effpid*effsel*effrec));
			d0Kinematics[1]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch/(effpid*effsel*effrec*effacc));

			//truth match to real, accepted and reconstructed D0
			if(D0TRUEIDX->at(s)<0) continue;
			int d = D0TRUEIDX->at(s);
			if(TRUEDTRK0IDX->at(d)==-1) continue;
			if(TRUEDTRK2IDX->at(d)!=-1) continue;
			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;
			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;
			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			fndD0recoPT.at(4)->Fill(JetPT,weight);
			fndD0recoPT.at(3)->Fill(JetPT,weight/effpid);
			fndD0recoPT.at(2)->Fill(JetPT,weight/(effpid*effsel));
			fndD0recoPT.at(1)->Fill(JetPT,weight/(effpid*effsel*effrec));
			fndD0recoPT.at(0)->Fill(JetPT,weight/(effpid*effsel*effrec*effacc));

			fndD0truePT.at(4)->Fill(JetTruePT,weight);
			fndD0truePT.at(3)->Fill(JetTruePT,weight/effpid);
			fndD0truePT.at(2)->Fill(JetTruePT,weight/(effpid*effsel));
			fndD0truePT.at(1)->Fill(JetTruePT,weight/(effpid*effsel*effrec));
			fndD0truePT.at(0)->Fill(JetTruePT,weight/(effpid*effsel*effrec*effacc));

			fndD02DPT.at(4)->Fill(JetTruePT,JetPT,weight);
			fndD02DPT.at(3)->Fill(JetTruePT,JetPT,weight/effpid);
			fndD02DPT.at(2)->Fill(JetTruePT,JetPT,weight/(effpid*effsel));
			fndD02DPT.at(1)->Fill(JetTruePT,JetPT,weight/(effpid*effsel*effrec));
			fndD02DPT.at(0)->Fill(JetTruePT,JetPT,weight/(effpid*effsel*effrec*effacc));

			foundD0s.push_back(d*100+s);

		//TODO	break;//only keep one D0 candidate per entry
		}
		for(unsigned int i=0; i<trueD0s.size(); ++i) {
			bool matched(false);
			for(unsigned int j=0; j<foundD0s.size(); ++j) {
				if(foundD0s[j]==trueD0s[i]) {
					matched=true;
					break;
				}
			}
			if(!matched) std::cout << progress.count()-1 << ": truD0 " << trueD0s[i] << " has no matching fndD0" << std::endl;
		}
		for(unsigned int i=0; i<foundD0s.size(); ++i) {
			bool matched(false);
			for(unsigned int j=0; j<trueD0s.size(); ++j) {
				if(foundD0s[i]==trueD0s[j]) {
					matched=true;
					break;
				}
			}
			if(!matched) std::cout << progress.count()-1 << ": fndD0 " << foundD0s[i] << " has no matching truD0" << std::endl;
		}
	}

	printf("bins of true jet pT\n");
	for(int i=0; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int j=1; j<=truD0truePT.at(i)->GetNbinsX(); ++j) {
			printf("%6.0f-%6.0f\t",truD0truePT.at(i)->GetBinLowEdge(j),truD0truePT.at(i)->GetBinLowEdge(j+1));
			printf("%5.0f\t%5.0f\t%5.0f\t%5.0f",truD0truePT.at(i)->GetBinContent(j),effD0truePT.at(i)->GetBinContent(j),fndD0truePT.at(i)->GetBinContent(j),fitD0truePT.at(i)->GetBinContent(j));
			if(i>0) {
				printf("\t%5.3f", effD0truePT.at(i-1)->GetBinContent(j)/truD0truePT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f\t%5.3f", truD0truePT.at(i)->GetBinContent(j)/truD0truePT.at(i-1)->GetBinContent(j),
						         fndD0truePT.at(i)->GetBinContent(j)/fndD0truePT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f", (fndD0truePT.at(i)->GetBinContent(j)/fndD0truePT.at(i-1)->GetBinContent(j))/
						  (truD0truePT.at(i)->GetBinContent(j)/truD0truePT.at(i-1)->GetBinContent(j)));
			}
			printf("\n");
		}
	}

	printf("bins of reco jet pT\n");
	for(int i=0; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int j=1; j<=truD0recoPT.at(i)->GetNbinsX(); ++j) {
			printf("%6.0f-%6.0f\t",truD0recoPT.at(i)->GetBinLowEdge(j),truD0recoPT.at(i)->GetBinLowEdge(j+1));
			printf("%5.0f\t%5.0f\t%5.0f\t%5.0f\t",truD0recoPT.at(i)->GetBinContent(j),effD0recoPT.at(i)->GetBinContent(j),fndD0recoPT.at(i)->GetBinContent(j),fitD0recoPT.at(i)->GetBinContent(j));
			if(i>0) {
				printf("\t%5.3f", effD0recoPT.at(i-1)->GetBinContent(j)/truD0recoPT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f\t%5.3f", truD0recoPT.at(i)->GetBinContent(j)/truD0recoPT.at(i-1)->GetBinContent(j),
						         fndD0recoPT.at(i)->GetBinContent(j)/fndD0recoPT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f", (fndD0recoPT.at(i)->GetBinContent(j)/fndD0recoPT.at(i-1)->GetBinContent(j))/
						  (truD0recoPT.at(i)->GetBinContent(j)/truD0recoPT.at(i-1)->GetBinContent(j)));
			}
			printf("\n");
		}
	}

	printf("bins of (X) true and (Y) reco jet pT\neff(found)/eff(truth)\n");
	for(int i=1; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int k=1; k<=truD02DPT.at(i)->GetNbinsX(); ++k) {
			printf("\t%3.0f-%3.0f\t",truD02DPT.at(i)->GetXaxis()->GetBinLowEdge(k)/1000.,truD02DPT.at(i)->GetXaxis()->GetBinLowEdge(k+1)/1000.);
		}
		printf("\n");
		for(int j=1; j<=truD02DPT.at(i)->GetNbinsY(); ++j) {
			printf("%3.0f-%3.0f",truD02DPT.at(i)->GetYaxis()->GetBinLowEdge(j)/1000.,truD02DPT.at(i)->GetYaxis()->GetBinLowEdge(j+1)/1000.);
			for(int k=1; k<=truD02DPT.at(i)->GetNbinsX(); ++k) {
				double ratio = (fndD02DPT.at(i)->GetBinContent(k,j)/fndD02DPT.at(i-1)->GetBinContent(k,j))/
				               (truD02DPT.at(i)->GetBinContent(k,j)/truD02DPT.at(i-1)->GetBinContent(k,j));
				double error = ratio*TMath::Sqrt(TMath::Power(fndD02DPT.at(i  )->GetBinError(k,j)/fndD02DPT.at(i  )->GetBinContent(k,j),2.)
						             //+ TMath::Power(fndD02DPT.at(i-1)->GetBinError(k,j)/fndD02DPT.at(i-1)->GetBinContent(k,j),2.)
						               + TMath::Power(truD02DPT.at(i  )->GetBinError(k,j)/truD02DPT.at(i  )->GetBinContent(k,j),2.));
						             //+ TMath::Power(truD02DPT.at(i-1)->GetBinError(k,j)/truD02DPT.at(i-1)->GetBinContent(k,j),2.));
				printf("\t%5.3f+/-%5.3f", ratio, error);
			}
			printf("\n");
		}
	}

	printf("bins of (X) true and (Y) reco jet pT\nN(found)\n");
	for(int i=1; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int k=1; k<=truD02DPT.at(i)->GetNbinsX(); ++k) {
			printf("\t%3.0f-%3.0f\t",truD02DPT.at(i)->GetXaxis()->GetBinLowEdge(k)/1000.,truD02DPT.at(i)->GetXaxis()->GetBinLowEdge(k+1)/1000.);
		}
		printf("\n");
		for(int j=1; j<=truD02DPT.at(i)->GetNbinsY(); ++j) {
			printf("%3.0f-%3.0f",truD02DPT.at(i)->GetYaxis()->GetBinLowEdge(j)/1000.,truD02DPT.at(i)->GetYaxis()->GetBinLowEdge(j+1)/1000.);
			for(int k=1; k<=truD02DPT.at(i)->GetNbinsX(); ++k) {
				printf("\t%6.0f\t", fndD02DPT.at(i)->GetBinContent(k,j));
			}
			printf("\n");
		}
	}

	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 0.00, 0.00};
	Double_t greens[NRGBs] = { 0.00, 0.50, 1.00, 0.50, 0.00};
	Double_t blues[NRGBs]  = { 0.00, 0.00, 1.00, 1.00, 1.00};
	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);
	TCanvas c1;
	c1.SetLogx();
	d0Kinematics[9]->Divide(d0Kinematics[7]);
	d0Kinematics[8]->Divide(d0Kinematics[6]);
	d0Kinematics[9]->Divide(d0Kinematics[8]);
	d0Kinematics[9]->SetMinimum(0.);
	d0Kinematics[9]->SetMaximum(2.);
	d0Kinematics[9]->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0KinePIDRatio.pdf");
	d0Kinematics[7]->Divide(d0Kinematics[5]);
	d0Kinematics[6]->Divide(d0Kinematics[4]);
	d0Kinematics[7]->Divide(d0Kinematics[6]);
	d0Kinematics[7]->SetMinimum(0.);
	d0Kinematics[7]->SetMaximum(2.);
	d0Kinematics[7]->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0KineSELRatio.pdf");
	d0Kinematics[5]->Divide(d0Kinematics[3]);
	d0Kinematics[4]->Divide(d0Kinematics[2]);
	d0Kinematics[5]->Divide(d0Kinematics[4]);
	d0Kinematics[5]->SetMinimum(0.);
	d0Kinematics[5]->SetMaximum(2.);
	d0Kinematics[5]->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0KineRECRatio.pdf");
	d0Kinematics[3]->Divide(d0Kinematics[1]);
	d0Kinematics[2]->Divide(d0Kinematics[0]);
	d0Kinematics[3]->Divide(d0Kinematics[2]);
	d0Kinematics[3]->SetMinimum(0.);
	d0Kinematics[3]->SetMaximum(2.);
	d0Kinematics[3]->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0KineACCRatio.pdf");
	//for(unsigned int i=0; i<d0Kinematics.size(); ++i) {
	//	d0Kinematics[i]->Draw("colz");
	//	TString plotname("D0Kinematics");
	//	plotname+=i;
	//	plotname+=".pdf";
	//	c1.SaveAs(gSaveDir+"/"+plotname);
	//}
	
	for(int i=1; i<=d0Kinematics[3]->GetNbinsX(); ++i) {
		for(int j=1; j<=d0Kinematics[3]->GetNbinsY(); ++j) {
			if(d0Kinematics[3]->GetBinContent(i,j)>0) {
				d0KinematicsPulls[0]->Fill((d0Kinematics[3]->GetBinContent(i,j)-1.)/d0Kinematics[3]->GetBinError(i,j));
			}
			if(d0Kinematics[5]->GetBinContent(i,j)>0) {
				d0KinematicsPulls[1]->Fill((d0Kinematics[5]->GetBinContent(i,j)-1.)/d0Kinematics[5]->GetBinError(i,j));
			}
			if(d0Kinematics[7]->GetBinContent(i,j)>0) {
				d0KinematicsPulls[2]->Fill((d0Kinematics[7]->GetBinContent(i,j)-1.)/d0Kinematics[7]->GetBinError(i,j));
			}
			if(d0Kinematics[9]->GetBinContent(i,j)>0) {
				d0KinematicsPulls[3]->Fill((d0Kinematics[9]->GetBinContent(i,j)-1.)/d0Kinematics[9]->GetBinError(i,j));
			}
		}
	}
	c1.SetLogx(0);
	d0KinematicsPulls[0]->Draw();
	c1.SaveAs(gSaveDir+"/D0KineACCPulls.pdf");
	d0KinematicsPulls[1]->Draw();
	c1.SaveAs(gSaveDir+"/D0KineRECPulls.pdf");
	d0KinematicsPulls[2]->Draw();
	c1.SaveAs(gSaveDir+"/D0KineSELPulls.pdf");
	d0KinematicsPulls[3]->Draw();
	c1.SaveAs(gSaveDir+"/D0KinePIDPulls.pdf");

	const Int_t NRGBs2 = 5;
	const Int_t NCont2 = 255;
	Double_t stops2[NRGBs2]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds2[NRGBs2]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
	Double_t greens2[NRGBs2] = { 1.00, 0.95, 0.50, 0.00, 0.00};
	Double_t blues2[NRGBs2]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);

	c1.Clear();
	c1.Divide(4,5);
	for(int i=0; i<5; ++i) {
		double max = 1.1*TMath::Max(TMath::Max(truD02DPT.at(i)->GetMaximum(),effD02DPT.at(i)->GetMaximum()),TMath::Max(fndD02DPT.at(i)->GetMaximum(),fitD02DPT.at(i)->GetMaximum()));
		truD02DPT.at(i)->SetMinimum(0.);
		truD02DPT.at(i)->SetMaximum(max);
		fndD02DPT.at(i)->SetMinimum(0.);
		fndD02DPT.at(i)->SetMaximum(max);
		fitD02DPT.at(i)->SetMinimum(0.);
		fitD02DPT.at(i)->SetMaximum(max);
		effD02DPT.at(i)->SetMinimum(0.);
		effD02DPT.at(i)->SetMaximum(max);
		c1.cd(4*i+1);
		truD02DPT.at(i)->Draw("colz");
		c1.cd(4*i+2);
		fndD02DPT.at(i)->Draw("colz");
		c1.cd(4*i+3);
		fitD02DPT.at(i)->Draw("colz");
		c1.cd(4*i+4);
		effD02DPT.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/2DEffMaps.pdf");

	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);
	c1.Clear();
	c1.Divide(5,3);
	for(int i=0; i<5; ++i) {
		fndD02DPT.at(i)->Divide(truD02DPT.at(i));
		fitD02DPT.at(i)->Divide(truD02DPT.at(i));
		effD02DPT.at(i)->Divide(truD02DPT.at(i));
		double max = 2.0;
		fndD02DPT.at(i)->SetMinimum(0.);
		fndD02DPT.at(i)->SetMaximum(max);
		fitD02DPT.at(i)->SetMinimum(0.);
		fitD02DPT.at(i)->SetMaximum(max);
		effD02DPT.at(i)->SetMinimum(0.);
		effD02DPT.at(i)->SetMaximum(max);
		c1.cd(i+1);
		fndD02DPT.at(i)->Draw("colz");
		c1.cd(i+6);
		fitD02DPT.at(i)->Draw("colz");
		c1.cd(i+11);
		effD02DPT.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/2DEffMapsRatios.pdf");

	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);
	c1.Clear();
	c1.Divide(6,7);
	for(uint i=0; i<42; ++i) {
		c1.cd(i+1);
		d0KinematicsTRBinned.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/D0KinematicsByJPTBin.pdf");
	for(uint i=0; i<42; ++i) {
		c1.cd(i+1);
		d0KinematicsTRBinned_recD_tru.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/D0KinematicsByJPTBin_recD_truKine.pdf");
	for(uint i=0; i<42; ++i) {
		c1.cd(i+1);
		d0KinematicsTRBinned_recD_rec.at(i)->Draw("colz");
	}
	c1.SaveAs(gSaveDir+"/D0KinematicsByJPTBin_recD_recKine.pdf");

	dm->reset();

	return true;
}

bool SimDFitter::fitD(double& yield4, double& error4, double& yield5, double& error5, uint binPT, uint binY, uint effType) {
	if(!_inputsSet) return false;
	if(!_effsSet) {
		std::cout << "INFO in SimDFitter::fitD: D efficiencies not calculated yet" << std::endl;
		std::cout << "                          running now" << std::endl;
		if(!addEffs()) return false;
	}

	//check if we've changed bin (need to recreate D pT bins)
	if(_jptmin != _ptBins[binPT] || _jptmax != _ptBins[binPT+1] ||
	   _ymin != _yBins[binY] || _ymax != _yBins[binY+1]) {
		_jptmin = _ptBins[binPT];
		_jptmax = _ptBins[binPT+1];
		_ymin = _yBins[binY];
		_ymax = _yBins[binY+1];
		_initialised = false;
	}

	if(effType==4) {
		_flavour = 4;
		_useEffs = true;
	} else if(effType==5) {
		_flavour = 5;
		_useEffs = true;
	} else {
		_useEffs = false;
	}
	
	if(!doFits()) return false;
	std::cout << _yield4 << "\t" << _error4 << "\t" << _yield5 << "\t" << _error5 << std::endl;//TODO

	yield5 = _yield5;
	error5 = _error5;
	yield4 = _yield4;
	error4 = _error4;

	return true;
}

//bool SimDFitter::makeInput(int flav) {
//	TString fname = "/data/dijets/dijets_";
//	if(flav==0) {
//		fname+="2016.root";
//	} else {
//		fname+="sim";
//		fname+=flav;
//		fname+=".root";
//	}
//
//	TFile* f = TFile::Open(fname);
//	if(!f) return false;
//
//	TTree* t = dynamic_cast<TTree*>(f->Get("T"));
//	if(!t) return false;
//
//	std::vector<double> *vD0M = new std::vector<double>();
//	std::vector<double> *vD0IPCHI2 = new std::vector<double>();
//	std::vector<double> *vD0PT = new std::vector<double>();
//	std::vector<double> *vD0PX = new std::vector<double>();
//	std::vector<double> *vD0PY = new std::vector<double>();
//	std::vector<double> *vD0PZ = new std::vector<double>();
//	std::vector<double> *vD0E = new std::vector<double>();
//	std::vector<double> *vD0X = new std::vector<double>();
//	std::vector<double> *vD0Y = new std::vector<double>();
//	std::vector<double> *vD0Z = new std::vector<double>();
//	std::vector<double> *vD0KP = new std::vector<double>();
//	std::vector<double> *vD0KPT = new std::vector<double>();
//	std::vector<double> *vD0KPX = new std::vector<double>();
//	std::vector<double> *vD0KPY = new std::vector<double>();
//	std::vector<double> *vD0KPZ = new std::vector<double>();
//	std::vector<double> *vD0PIP = new std::vector<double>();
//	std::vector<double> *vD0PIPT = new std::vector<double>();
//	std::vector<double> *vD0PIPX = new std::vector<double>();
//	std::vector<double> *vD0PIPY = new std::vector<double>();
//	std::vector<double> *vD0PIPZ = new std::vector<double>();
//	std::vector<double> *vD0KPNNK = new std::vector<double>();
//	std::vector<double> *vD0PIPNNPI = new std::vector<double>();
//	std::vector<double> *vD0KWEIGHT = new std::vector<double>();
//	std::vector<double> *vD0PIWEIGHT = new std::vector<double>();
//	std::vector<double> *vD0TRUEIDX = new std::vector<double>();
//	std::vector<double> *vD0FROMB = new std::vector<double>();
//
//	double JetPT;
//	double JetEta;
//	double JetTruePT;
//
//	t->SetBranchAddress("D0M",           &vD0M);
//	t->SetBranchAddress("D0IPCHI2",      &vD0IPCHI2);
//	t->SetBranchAddress("D0PT",          &vD0PT);
//	t->SetBranchAddress("D0PX",          &vD0PX);
//	t->SetBranchAddress("D0PY",          &vD0PY);
//	t->SetBranchAddress("D0PZ",          &vD0PZ);
//	t->SetBranchAddress("D0E",           &vD0E);
//	t->SetBranchAddress("D0X",           &vD0X);
//	t->SetBranchAddress("D0Y",           &vD0Y);
//	t->SetBranchAddress("D0Z",           &vD0Z);
//	t->SetBranchAddress("D0KP",          &vD0KP);
//	t->SetBranchAddress("D0KPT",         &vD0KPT);
//	t->SetBranchAddress("D0KPX",         &vD0KPX);
//	t->SetBranchAddress("D0KPY",         &vD0KPY);
//	t->SetBranchAddress("D0KPZ",         &vD0KPZ);
//	t->SetBranchAddress("D0PIP",         &vD0PIP);
//	t->SetBranchAddress("D0PIPT",        &vD0PIPT);
//	t->SetBranchAddress("D0PIPX",        &vD0PIPX);
//	t->SetBranchAddress("D0PIPY",        &vD0PIPY);
//	t->SetBranchAddress("D0PIPZ",        &vD0PIPZ);
//	t->SetBranchAddress("D0KPNNK",       &vD0KPNNK);
//	t->SetBranchAddress("D0PIPNNPI",     &vD0PIPNNPI);
//	t->SetBranchAddress("D0KWEIGHT",     &vD0KWEIGHT);
//	t->SetBranchAddress("D0PIWEIGHT",    &vD0PIWEIGHT);
//	if(flav!=0) {
//		t->SetBranchAddress("D0TRUEIDX",     &vD0TRUEIDX);
//		t->SetBranchAddress("D0FROMB",       &vD0FROMB);
//	}
//
//	t->SetBranchAddress("JetPT",         &JetPT);
//	t->SetBranchAddress("JetEta",        &JetEta);
//	if(flav!=0) {
//		t->SetBranchAddress("JetTruePT",     &JetTruePT);
//	}
//
//	unsigned int nentries0 = t->GetEntries();
//
//	TFile* fout = TFile::Open("d0InputDatasets.root","UPDATE");
//	TString tname = "T";
//	tname+=flav;
//	TTree* tout = new TTree(tname,"");
//
//	double D0M(0.), D0PT(0.), D0Eta(0.), D0LogIPChi2(0.);
//
//	tout->Branch("JetPT",        &JetPT);
//	tout->Branch("JetEta",       &JetEta);
//	tout->Branch("JetTruePT",    &JetTruePT);
//	tout->Branch("D0M",          &D0M);
//	tout->Branch("D0PT",         &D0PT);
//	tout->Branch("D0Eta",        &D0Eta);
//	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
//
//	//first make non-vector tree for fits
//	boost::progress_display progress( nentries0 );
//	//for(unsigned int ientry=0; ientry<nentries0; ++ientry)
//	for(uint i=0; i<nentries0; ++i) {
//		t->GetEntry(i);
//		++progress;
//		//t0->GetEntry(ientry);
//
//		for(unsigned int s=0; s<vD0M->size(); ++s) {
//			//first check in our tight acceptance
//			TVector3 D0P (vD0PX->at(s)  ,vD0PY->at(s)  ,vD0PZ->at(s));
//			TVector3 D0P0(vD0KPX->at(s) ,vD0KPY->at(s) ,vD0KPZ->at(s));
//			TVector3 D0P1(vD0PIPX->at(s),vD0PIPY->at(s),vD0PIPZ->at(s));
//
//			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
//			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
//			//if(!(D0P.Pt()>_dptmin)) continue;
//			//if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;
//			if(flav!=0 && vD0TRUEIDX->at(s)==-1) continue;
//			if(flav==5 && vD0FROMB->at(s)==0) continue;
//			if(flav==4 && vD0FROMB->at(s)==1) continue;
//			if(flav==0 && vD0M->at(s)>1844. && vD0M->at(s)<1884.) continue;
//
//			//check PID cuts
//			//if(!_dataIsMC && vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
//			//if(!dataIsMC && vD0PIPNNPI->at(s)<0.1 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //pion PID turned off
//
//			D0M         = vD0M->at(s);
//			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
//			D0PT        = D0P.Pt();
//			D0Eta       = D0P.Eta();
//
//			//for the combinatorial, we need to resample the mass to remove signal window
//			if(flav==0) D0M = gRandom->Uniform(1784.,1944.);
//
//			tout->Fill();
//			break;//only keep one D0 candidate per entry
//		}
//	}
//
//	tout->AutoSave();
//	fout->Close();
//
//	return true;
//}

bool SimDFitter::loadDatasets() {
	TFile* f = TFile::Open(dFileName());
	//if(!f) {
	//	if(!makeInput(0)) return false;
	//	if(!makeInput(4)) return false;
	//	if(!makeInput(5)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//}
	if(!f) return false;

	TTree* t0 = dynamic_cast<TTree*>(f->Get("T0"));
	//if(!t0) {
	//	if(!makeInput(0)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//	t0 = dynamic_cast<TTree*>(f->Get("T0"));
	//}
	TTree* t4 = dynamic_cast<TTree*>(f->Get("T4"));
	//if(!t4) {
	//	if(!makeInput(4)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//	t4 = dynamic_cast<TTree*>(f->Get("T4"));
	//}
	TTree* t5 = dynamic_cast<TTree*>(f->Get("T5"));
	//if(!t5) {
	//	if(!makeInput(5)) return false;
	//	f = TFile::Open("d0InputDatasets.root");
	//	t5 = dynamic_cast<TTree*>(f->Get("T5"));
	//}
	TTree* t7 = dynamic_cast<TTree*>(f->Get("T7"));
	if(!t0 || !t4 || !t5 || !t7) return false;

	//setup equal bins in pT(D)
	setupDPtBins(t7);

	RooRealVar D0M        ("D0M",         "D0M",       1784., 1944., "MeV/#it{c}^{2}"); 
	RooRealVar D0LOGIPCHI2("D0LogIPChi2", "D0LogIPChi2", -5.,   15., ""); 
	RooRealVar D0PT       ("D0PT",        "D0PT",         0.,   1e5, "MeV/#it{c}"); 
	RooRealVar JetPT      ("JetPT",       "JetPT",        0.,   1e5, "MeV/#it{c}"); 
	RooRealVar JetEta     ("JetEta",      "JetEta",       0.,   10., ""); 
	RooRealVar weight4    ("weight4",      "weight4",     0.,   30., ""); 
	RooRealVar weight5    ("weight5",      "weight5",     0.,   30., ""); 

	RooThresholdCategory* ptCat(0);

	TString binName="PTbin";
	binName+=_ptBinsD.size()-2;
	if(!_usePtFracBins) {
		ptCat = new RooThresholdCategory("ptCat","ptCat", D0PT, binName, _ptBinsD.size()-2);
		for(uint i=1; i<_ptBinsD.size()-1; ++i) {
			binName="PTbin";
			binName+=i-1;
			ptCat->addThreshold(_ptBinsD[i],binName.Data(),i-1);
		}
	} else {
		RooFormulaVar* D0PTFrac = new RooFormulaVar("D0PTFrac","","@0/@1",RooArgList(D0PT,JetPT));
		ptCat = new RooThresholdCategory("ptCat","ptCat", *D0PTFrac, binName, _nDPtBins-1);
		for(uint i=1; i<_nDPtBins; ++i) {
			binName="PTbin";
			binName+=i-1;
			ptCat->addThreshold(static_cast<double>(i)/_nDPtBins,binName.Data(),i-1);
		}
	}

	binName="JetPTbin";
	binName+=_ptBins.size()-2;
	RooThresholdCategory jetptCat("jetptCat","jetptCat", JetPT, binName, _ptBins.size()-2);
	for(uint i=1; i<_ptBins.size()-1; ++i) {
		binName="JetPTbin";
		binName+=i-1;
		jetptCat.addThreshold(_ptBins[i],binName.Data(),i-1);
	}

	binName="JetYbin";
	binName+=_yBins.size()-2;
	RooThresholdCategory jetyCat("jetyCat","jetyCat", JetEta, binName, _yBins.size()-2);
	for(uint i=1; i<_yBins.size()-1; ++i) {
		binName="JetYbin";
		binName+=i-1;
		jetyCat.addThreshold(_yBins[i],binName.Data(),i-1);
	}

	RooCategory typeCat("typeCat","typeCat");
	typeCat.defineType("data",3);
	typeCat.defineType("charm",1);
	typeCat.defineType("beauty",2);
	typeCat.defineType("sideband",0);

	RooDataSet cData("cData", "cData", RooArgList(D0LOGIPCHI2,D0M,D0PT,JetPT,JetEta,weight4,weight5), RooFit::Import(*t4));
	RooDataSet bData("bData", "bData", RooArgList(D0LOGIPCHI2,D0M,D0PT,JetPT,JetEta,weight4,weight5), RooFit::Import(*t5));
	RooDataSet sData("sData", "sData", RooArgList(D0LOGIPCHI2,D0M,D0PT,JetPT,JetEta,weight4,weight5), RooFit::Import(*t0));
	RooDataSet dData("dData", "dData", RooArgList(D0LOGIPCHI2,D0M,D0PT,JetPT,JetEta,weight4,weight5), RooFit::Import(*t7));
	RooDataSet data("data", "data", RooArgList(D0LOGIPCHI2,D0M,D0PT,JetPT,JetEta,weight4,weight5), RooFit::Index(typeCat), RooFit::Import("charm",cData), RooFit::Import("beauty",bData), RooFit::Import("sideband",sData), RooFit::Import("data",dData));
	data.addColumn(*ptCat);
	data.addColumn(jetptCat);
	data.addColumn(jetyCat);
	sData.addColumn(*ptCat);

	ws->import(cData);
	ws->import(bData);
	ws->import(sData);
	ws->import(dData);
	ws->import(data);
	nEntries = data.sumEntries();

	return true;
}

bool SimDFitter::setupDPtBins(TTree* t) {
	//TODO
	//_ptBinsD.clear();
	//_ptBinsD.push_back(0);
	//_ptBinsD.push_back(7e3);
	//_ptBinsD.push_back(10e3);
	//_ptBinsD.push_back(15e3);
	//_ptBinsD.push_back(20e3);
	//_ptBinsD.push_back(100e3);
	//return true;
	//TODO
	uint nbins(500);
	TH1D pt_pdf("pt_pdf","",nbins,0.,100e3);
	//TH1D pt_cdf("pt_cdf","",nbins,0.,100e3);

	TString cutStr = TString::Format("JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f",_jptmin,_jptmax,_ymin,_ymax);
	t->Draw("D0PT>>pt_pdf",cutStr);
	double total = pt_pdf.GetEntries();
	double sum(0.);
	double step = total/_nDPtBins;

	_ptBinsD.clear();
	_ptBinsD.push_back(0.);
	for(uint iBin=1; iBin<=nbins; ++iBin) {
		sum += pt_pdf.GetBinContent(iBin);
		if(sum >= _ptBinsD.size()*step) {
			std::cout << iBin << " " << pt_pdf.GetBinLowEdge(iBin) << " " << sum << " " << total << " " << step << " " << _ptBinsD.size()*step << std::endl;//TODO
			_ptBinsD.push_back(pt_pdf.GetBinLowEdge(iBin+1));
		}
		//pt_cdf.SetBinContent(i,sum);
	}
	//pt_cdf.Scale(1./sum);
	_ptBinsD.back() = 1e5;//TODO

	//TCanvas c;
	//pt_cdf.Draw();
	//c.SaveAs("test.pdf");//TODO
	for(uint iBin=0; iBin<_ptBinsD.size(); ++iBin) {
		std::cout << _ptBinsD[iBin] << std::endl;
	}

	return true;
}

bool SimDFitter::setupModel() {
	RooRealVar* DM = ws->var("D0M");
	RooRealVar* DLOGIPCHI2 = ws->var("D0LogIPChi2");

	//mass part
	RooRealVar dMean("dMean","$\\mu_m$",valMassMean,1844.,1884.);
	RooRealVar dBaseWidth("dBaseWidth","$\\sigma_C$",valMassWidth,4.,40.);
	RooRealVar dWidthScale("dWidthScale","$\\sigma_C^{\\rm data}/\\sigma_C^{\\rm sim}$",1.,0.5,2.);
	RooFormulaVar dWidth("dWidth","","@0*@1",RooArgList(dBaseWidth,dWidthScale));

	RooRealVar dRatio("dRatio","$r$", valMassRatio, 1.0, 5);
	RooFormulaVar dWidthG("dWidthG","","@0*@1",RooArgList(dWidth,dRatio));
	RooRealVar dAlpha("dAlpha", "$\\alpha$", valMassAlpha, 0.5, 5.0);
	RooRealVar dN(    "dN",     "$N$",     valMassN);
	RooRealVar dFrac("dFrac","$f$", valMassFrac, 0., 1.);

	RooGaussian signalMass_gauss("signalMass_gauss","",*DM,dMean,dWidthG);
	RooCBShape  signalMass_cb("signalMass_cb", "", *DM, dMean, dWidth, dAlpha, dN); 
	RooAddPdf signalMass( "signalMass", "", RooArgList(signalMass_cb,signalMass_gauss), RooArgList(dFrac) );

	RooRealVar p0("p0","$p_0$",0., -1., 10.);
	RooPolynomial bkgrndMass("bkgrndMass","",*DM, RooArgList(p0));

	//IP part
	RooRealVar promptBaseMean("promptBaseMean",  "$\\mu^{\\rm prompt}$",     valPromptMean,0.5,1.);
	RooRealVar promptMeanShift("promptMeanShift","$\\mu^{\\rm prompt data}/\\mu^{\\rm prompt sim}$",0.,-1.,1.);
	RooFormulaVar promptMean("promptMean","","@0+@1",RooArgList(promptBaseMean,promptMeanShift));
	RooRealVar promptWidth("promptWidth",        "$\\sigma^{\\rm prompt}$",  valPromptWidth,0.5,1.5);
	RooRealVar promptAsym("promptAsym",          "$\\epsilon^{\\rm prompt}$",valPromptAsym,-0.5,0.);
	RooRealVar promptRhoL("promptRhoL",          "$\\rho_L^{\\rm prompt}$",  valPromptRhoL,0.5,2.);
	RooRealVar promptRhoR("promptRhoR",          "$\\rho_R^{\\rm prompt}$",  valPromptRhoR);//,0.5,2.);
	RooPromptShape promptIP("promptIP","",*DLOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	RooRealVar displcMean("displcMean",  "$\\mu^{\\rm displc}$",     valDisplMean,2.,12.);
	RooRealVar displcWidth("displcWidth","$\\sigma^{\\rm displc}$",  valDisplWidth,0.5,6.);
	RooRealVar displcAsym("displcAsym",  "$\\epsilon^{\\rm displc}$",valDisplAsym,-1.,1.);
	RooRealVar displcRhoL("displcRhoL",  "$\\rho_L^{\\rm displc}$",  valDisplRhoL);
	RooRealVar displcRhoR("displcRhoR",  "$\\rho_R^{\\rm displc}$",  valDisplRhoR);
	RooPromptShape displcIP("displcIP","",*DLOGIPCHI2,displcMean,displcWidth,displcAsym,displcRhoL,displcRhoR);

	RooKeysPdf bkgrndIP("bkgrndIP","",*DLOGIPCHI2,*dynamic_cast<RooDataSet*>(ws->data("sData")));
//TODO//20200420	RooDataSet* sData = dynamic_cast<RooDataSet*>(ws->data("sData"));
//TODO//20200420	RooRealVar bkgrndBinIdx("bkgrndBinIdx","pt bin index",0.,-1.,_ptBinsD.size());
//TODO//20200420	RooArgList pdfs("pdfs");
//TODO//20200420	const RooArgSet* vars = sData->get();
//TODO//20200420	for(uint iBin=0; iBin<_ptBinsD.size()-1; ++iBin) {
//TODO//20200420		TString nameStr = TString::Format("_%d",iBin);
//TODO//20200420		TString cutStr = TString::Format("ptCat==%d",iBin);
//TODO//20200420		RooDataSet* ds = new RooDataSet("sData"+nameStr, "", sData, *vars, cutStr);
//TODO//20200420		RooKeysPdf* pdf = new RooKeysPdf("bkgrndIP"+nameStr, "", *DLOGIPCHI2, *ds);
//TODO//20200420		pdfs.add(*pdf);
//TODO//20200420	}
//TODO//20200420	RooMultiKeysPdf bkgrndIP("bkgrndIP", "", *DLOGIPCHI2, bkgrndBinIdx, pdfs);

	RooProdPdf prompt("prompt","",RooArgSet(signalMass,promptIP));
	RooProdPdf displc ("displc", "",RooArgSet(signalMass,displcIP));
	RooProdPdf bkgrnd("bkgrnd","",RooArgSet(bkgrndMass,bkgrndIP));

	// -- yields
	RooRealVar signalYield( "signalYield", "$N_{\\rm sig}$",   -10.,  nEntries);
	RooRealVar promptYield( "promptYield", "$N_{\\rm prmpt}$", -10.,  nEntries);
	RooRealVar displcYield( "displcYield", "$N_{\\rm displ}$", -10.,  nEntries);
	RooRealVar bkgrndYield( "bkgrndYield", "$N_{\\rm bkgnd}$", -10.,  nEntries);

	RooAddPdf mass_pdf( "mass_pdf",  "mass_pdf", RooArgList(signalMass,bkgrndMass), RooArgList(signalYield,bkgrndYield) );
	RooAddPdf ip_pdf(   "ip_pdf",    "ip_pdf",   RooArgList(promptIP,displcIP,bkgrndIP), RooArgList(promptYield,displcYield,bkgrndYield) );
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(prompt,displc,bkgrnd), RooArgList(promptYield,displcYield,bkgrndYield) );

	ws->importClassCode(RooPromptShape::Class());
	ws->importClassCode(RooKeysPdf::Class());
	ws->importClassCode(RooMultiKeysPdf::Class());
	ws->import(mass_pdf);
	ws->import(ip_pdf);
	ws->import(data_pdf,RooFit::RenameConflictNodes("_2d"));//RooFit::RecycleConflictNodes());

	RooSimWSTool sct(*ws);
	sct.build("sim_mass_pdf","mass_pdf",
	  	  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,signalYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	  	  RooFit::SplitParam("dBaseWidth,dRatio,dAlpha,dFrac,signalYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	sct.build("sim_ip_pdf","ip_pdf",
//TODO//20200420	  	  RooFit::SplitParam("promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,bkgrndBinIdx,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	  	  RooFit::SplitParam("promptBaseMean,promptWidth,promptAsym,promptRhoL,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
	sct.build("sim_data_pdf","data_pdf",
		  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
//TODO//20200420		  RooFit::SplitParam("dWidth,dRatio,dAlpha,dFrac,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,bkgrndBinIdx,promptYield,displcYield,bkgrndYield,p0,dMean","ptCat"));//,jetptCat,jetyCat"));
		  RooFit::SplitParam("dBaseWidth,dRatio,dAlpha,dFrac,promptBaseMean,promptWidth,promptAsym,promptRhoL,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));

//TODO//20200420	for(uint iBin=0; iBin<_ptBinsD.size()-1; ++iBin) {
//TODO//20200420		RooRealVar* binIdx = ws->var(TString::Format("bkgrndBinIdx_PTbin%d",iBin));
//TODO//20200420		binIdx->setConstant(true);
//TODO//20200420		binIdx->setVal(iBin);
//TODO//20200420	}

	//ws->Print();
	//TODO//initPars();
	//return false;

	return true;
}

//bool SimDFitter::setupDPtBins() {
//	
//	TString cutStr = TString::Format("JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f",_jptmin,_jptmax,_ymin,_ymax);
//	TString nameStr = TString::Format("_%.1f-%.1f_%.3f-%.3f",_jptmin,_jptmax,_ymin,_ymax);
//
//	RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data"));
//	RooRealVar* D0PT = ws->var("D0PT");
//	const RooArgSet* vars = data->get();
//	RooDataSet* ds = new RooDataSet("data"+nameStr, "", data, *vars, cutStr);
//
//	TString binName="PTbin";
//	binName+=_ptBinsD.size()-2;
//	RooThresholdCategory ptCat("ptCat"+nameStr,"ptCat", *D0PT, binName, _ptBinsD.size()-2);
//	for(uint i=1; i<_ptBinsD.size()-1; ++i) {
//		binName="PTbin";
//		binName+=i-1;
//		ptCat.addThreshold(_ptBinsD[i],binName.Data(),i-1);
//	}
//
//	ds->addColumn(ptCat);
//	RooSimWSTool sct(*ws);
//	sct.build("sim_mass_pdf"+nameStr,"mass_pdf",
//	  	  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,signalYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
//	  	  RooFit::SplitParam("dWidth,dRatio,dAlpha,dFrac,signalYield,bkgrndYield","ptCat"+nameStr));
//	sct.build("sim_ip_pdf"+nameStr,"ip_pdf",
//	  	  RooFit::SplitParam("promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"+nameStr));
//	sct.build("sim_data_pdf"+nameStr,"data_pdf",
//		  //RooFit::SplitParam("dMean,dWidth,dRatio,dAlpha,dFrac,p0,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"));//,jetptCat,jetyCat"));
//		  RooFit::SplitParam("dWidth,dRatio,dAlpha,dFrac,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR,displcMean,displcWidth,displcAsym,promptYield,displcYield,bkgrndYield","ptCat"+nameStr));
//
//	ws->Print();
//	//TODO//initPars();
//
//	return true;
//}

void SimDFitter::initPars() {
//signal mass
	ws->var(         "dAlpha_PTbin0"  )->setVal(  1.9180e+00);// +/-  6.72e-02
	ws->var(         "dAlpha_PTbin1"  )->setVal(  1.8576e+00);// +/-  4.03e-02
	ws->var(         "dAlpha_PTbin2"  )->setVal(  2.0967e+00);// +/-  6.33e-02
	ws->var(         "dAlpha_PTbin3"  )->setVal(  2.0526e+00);// +/-  8.37e-02
	ws->var(         "dAlpha_PTbin4"  )->setVal(  2.1923e+00);// +/-  1.05e-01
	ws->var(          "dFrac_PTbin0"  )->setVal(  9.3494e-01);// +/-  1.24e-02
	ws->var(          "dFrac_PTbin1"  )->setVal(  9.4841e-01);// +/-  8.24e-03
	ws->var(          "dFrac_PTbin2"  )->setVal(  9.2646e-01);// +/-  8.37e-03
	ws->var(          "dFrac_PTbin3"  )->setVal(  8.8369e-01);// +/-  1.61e-02
	ws->var(          "dFrac_PTbin4"  )->setVal(  7.9370e-01);// +/-  1.95e-02
	ws->var(          "dMean_PTbin0"  )->setVal(  1.8656e+03);// +/-  1.44e-01
	ws->var(          "dMean_PTbin1"  )->setVal(  1.8655e+03);// +/-  1.01e-01
	ws->var(          "dMean_PTbin2"  )->setVal(  1.8655e+03);// +/-  1.03e-01
	ws->var(          "dMean_PTbin3"  )->setVal(  1.8654e+03);// +/-  1.47e-01
	ws->var(          "dMean_PTbin4"  )->setVal(  1.8655e+03);// +/-  1.43e-01
	ws->var(         "dRatio_PTbin0"  )->setVal(  3.5897e+00);// +/-  3.58e-01
	ws->var(         "dRatio_PTbin1"  )->setVal(  3.4258e+00);// +/-  3.14e-01
	ws->var(         "dRatio_PTbin2"  )->setVal(  3.5701e+00);// +/-  2.29e-01
	ws->var(         "dRatio_PTbin3"  )->setVal(  2.8388e+00);// +/-  1.73e-01
	ws->var(         "dRatio_PTbin4"  )->setVal(  2.4915e+00);// +/-  8.27e-02
	ws->var(         "dWidth_PTbin0"  )->setVal(  7.5904e+00);// +/-  1.47e-01
	ws->var(         "dWidth_PTbin1"  )->setVal(  8.0534e+00);// +/-  1.09e-01
	ws->var(         "dWidth_PTbin2"  )->setVal(  8.9742e+00);// +/-  1.11e-01
	ws->var(         "dWidth_PTbin3"  )->setVal(  9.5189e+00);// +/-  1.69e-01
	ws->var(         "dWidth_PTbin4"  )->setVal(  1.1444e+01);// +/-  2.16e-01
	//ws->var(    "signalYield_PTbin0"  )->setVal(  1.4054e+04);// +/-  2.28e+02
	//ws->var(    "signalYield_PTbin1"  )->setVal(  2.7708e+04);// +/-  3.00e+02
	//ws->var(    "signalYield_PTbin2"  )->setVal(  3.2729e+04);// +/-  3.19e+02
	//ws->var(    "signalYield_PTbin3"  )->setVal(  1.8350e+04);// +/-  2.31e+02
	//ws->var(    "signalYield_PTbin4"  )->setVal(  3.6071e+04);// +/-  3.39e+02
//bkgrnd mass
	//ws->var(    "bkgrndYield_PTbin0"  )->setVal(  4.8410e+03);// +/-  1.31e+02
	//ws->var(    "bkgrndYield_PTbin1"  )->setVal(  7.6104e+03);// +/-  1.55e+02
	//ws->var(    "bkgrndYield_PTbin2"  )->setVal(  7.6209e+03);// +/-  1.50e+02
	//ws->var(    "bkgrndYield_PTbin3"  )->setVal(  4.0486e+03);// +/-  1.09e+02
	//ws->var(    "bkgrndYield_PTbin4"  )->setVal(  5.3398e+03);// +/-  1.29e+02
	ws->var(             "p0_PTbin0"  )->setVal( -4.8470e-04);// +/-  4.24e-06
	ws->var(             "p0_PTbin1"  )->setVal( -4.7094e-04);// +/-  5.41e-06
	ws->var(             "p0_PTbin2"  )->setVal( -4.6175e-04);// +/-  6.86e-06
	ws->var(             "p0_PTbin3"  )->setVal( -4.5688e-04);// +/-  1.07e-05
	ws->var(             "p0_PTbin4"  )->setVal( -4.1021e-04);// +/-  2.87e-05
//prompt IP
	ws->var(     "promptAsym_PTbin0"  )->setVal( -2.44938e-01);// +/-  1.94264e-02
	ws->var(     "promptAsym_PTbin1"  )->setVal( -3.85341e-01);// +/-  1.12410e-02
	ws->var(     "promptAsym_PTbin2"  )->setVal( -2.32890e-01);// +/-  1.20989e-02
	ws->var(     "promptAsym_PTbin3"  )->setVal( -1.58127e-01);// +/-  1.57671e-02
	ws->var(     "promptAsym_PTbin4"  )->setVal( -2.20983e-01);// +/-  1.23541e-02
	ws->var(     "promptMean_PTbin0"  )->setVal(  7.41476e-01);// +/-  3.22027e-02
	ws->var(     "promptMean_PTbin1"  )->setVal(  8.69684e-01);// +/-  2.10475e-02
	ws->var(     "promptMean_PTbin2"  )->setVal(  6.84181e-01);// +/-  1.58891e-02
	ws->var(     "promptMean_PTbin3"  )->setVal(  6.49780e-01);// +/-  2.59778e-02
	ws->var(     "promptMean_PTbin4"  )->setVal(  7.91027e-01);// +/-  1.74429e-02
	ws->var(     "promptRhoL_PTbin0"  )->setVal(  1.26607e+00);// +/-  4.62087e-02
	ws->var(     "promptRhoL_PTbin1"  )->setVal(  1.63344e+00);// +/-  5.11881e-02
	ws->var(     "promptRhoL_PTbin2"  )->setVal(  1.23796e+00);// +/-  3.79378e-02
	ws->var(     "promptRhoL_PTbin3"  )->setVal(  1.12382e+00);// +/-  2.89408e-02
	ws->var(     "promptRhoL_PTbin4"  )->setVal(  1.32783e+00);// +/-  3.59866e-02
	ws->var(     "promptRhoR_PTbin0"  )->setVal(  1.90330e+00);// +/-  5.90527e-02
	ws->var(     "promptRhoR_PTbin1"  )->setVal(  1.56755e+00);// +/-  3.65293e-02
	ws->var(     "promptRhoR_PTbin2"  )->setVal(  1.65797e+00);// +/-  3.21594e-02
	ws->var(     "promptRhoR_PTbin3"  )->setVal(  1.87261e+00);// +/-  5.23660e-02
	ws->var(     "promptRhoR_PTbin4"  )->setVal(  1.73098e+00);// +/-  4.19735e-02
	ws->var(    "promptWidth_PTbin0"  )->setVal(  1.12137e+00);// +/-  1.50698e-02
	ws->var(    "promptWidth_PTbin1"  )->setVal(  1.17254e+00);// +/-  8.66453e-03
	ws->var(    "promptWidth_PTbin2"  )->setVal(  1.09664e+00);// +/-  1.11622e-02
	ws->var(    "promptWidth_PTbin3"  )->setVal(  1.07137e+00);// +/-  1.07803e-02
	ws->var(    "promptWidth_PTbin4"  )->setVal(  1.15567e+00);// +/-  9.26588e-03
	//ws->var(    "promptYield_PTbin0"  )->setVal(  2.31630e+03);// +/-  2.69955e+01
	//ws->var(    "promptYield_PTbin1"  )->setVal(  8.68934e+03);// +/-  5.83562e+01
	//ws->var(    "promptYield_PTbin2"  )->setVal(  1.25728e+04);// +/-  7.68802e+01
	//ws->var(    "promptYield_PTbin3"  )->setVal(  7.29075e+03);// +/-  5.99772e+01
	//ws->var(    "promptYield_PTbin4"  )->setVal(  1.87655e+04);// +/-  8.82830e+01
//disp IP
	ws->var(     "displcAsym_PTbin0"  )->setVal(  1.1835e-02);// +/-  3.07e-02
	ws->var(     "displcAsym_PTbin1"  )->setVal( -7.3169e-02);// +/-  2.52e-02
	ws->var(     "displcAsym_PTbin2"  )->setVal( -2.1166e-01);// +/-  2.45e-02
	ws->var(     "displcAsym_PTbin3"  )->setVal( -1.6608e-01);// +/-  3.28e-02
	ws->var(     "displcAsym_PTbin4"  )->setVal( -2.5108e-01);// +/-  2.66e-02
	ws->var(     "displcMean_PTbin0"  )->setVal(  3.6616e+00);// +/-  1.23e-01
	ws->var(     "displcMean_PTbin1"  )->setVal(  4.1797e+00);// +/-  1.02e-01
	ws->var(     "displcMean_PTbin2"  )->setVal(  4.8858e+00);// +/-  9.70e-02
	ws->var(     "displcMean_PTbin3"  )->setVal(  4.7789e+00);// +/-  1.29e-01
	ws->var(     "displcMean_PTbin4"  )->setVal(  4.8852e+00);// +/-  1.01e-01
	ws->var(    "displcWidth_PTbin0"  )->setVal(  2.2864e+00);// +/-  2.94e-02
	ws->var(    "displcWidth_PTbin1"  )->setVal(  2.3217e+00);// +/-  2.31e-02
	ws->var(    "displcWidth_PTbin2"  )->setVal(  2.3108e+00);// +/-  2.44e-02
	ws->var(    "displcWidth_PTbin3"  )->setVal(  2.2657e+00);// +/-  2.97e-02
	ws->var(    "displcWidth_PTbin4"  )->setVal(  2.2072e+00);// +/-  2.47e-02
	//ws->var(    "displcYield_PTbin0"  )->setVal(  1.1203e+04);// +/-  2.04e+02
	//ws->var(    "displcYield_PTbin1"  )->setVal(  1.6159e+04);// +/-  2.28e+02
	//ws->var(    "displcYield_PTbin2"  )->setVal(  1.5315e+04);// +/-  2.20e+02
	//ws->var(    "displcYield_PTbin3"  )->setVal(  8.4982e+03);// +/-  1.57e+02
	//ws->var(    "displcYield_PTbin4"  )->setVal(  1.2530e+04);// +/-  1.97e+02
//bkgrnd IP
	//ws->var(    "bkgrndYield_PTbin0"  )->setVal(  4.8410e+03);// +/-  1.31e+02
	//ws->var(    "bkgrndYield_PTbin1"  )->setVal(  7.6104e+03);// +/-  1.55e+02
	//ws->var(    "bkgrndYield_PTbin2"  )->setVal(  7.6209e+03);// +/-  1.50e+02
	//ws->var(    "bkgrndYield_PTbin3"  )->setVal(  4.0486e+03);// +/-  1.09e+02
	//ws->var(    "bkgrndYield_PTbin4"  )->setVal(  5.3398e+03);// +/-  1.29e+02
}

bool SimDFitter::fit(fitType whichFit, fitSet whichSet, uint ptBinMin, uint ptBinMax, unsigned int nToys) {
	bool status=true;

	RooRealVar* DM = ws->var("D0M");
	RooRealVar* DIP = ws->var("D0LogIPChi2");

	TString name="";
	TString typeName="";
	TString cutStr="ptCat>=";
	cutStr+=ptBinMin;
	cutStr+=" && ptCat<";
	cutStr+=ptBinMax;
	cutStr+=" && JetPT>=";
	cutStr+=_jptmin;
	cutStr+=" && JetPT<";
	cutStr+=_jptmax;

	TString weightName="";
	if(_useEffs) {
		if(_flavour==4) weightName="weight4";
		if(_flavour==5) weightName="weight5";
	}

	switch(whichSet) {
		case fitSet::fitJustCharm:
			cutStr+=" && typeCat==1";
			name="_prompt";
			break;
		case fitSet::fitJustBeauty:
			cutStr+=" && typeCat==2";
			name="_displ";
			break;
		case fitSet::fitJustSignal:
			cutStr+=" && typeCat>0 && typeCat<3";
			name="_signal";
			break;
		case fitSet::fitJustSideband:
			cutStr+=" && typeCat==0";
			name="_bkgrnd";
			break;
		case fitSet::fitData:
			cutStr+=" && JetEta>=";
			cutStr+=_ymin;
			cutStr+=" && JetEta<";
			cutStr+=_ymax;
			cutStr+=" && typeCat==3";
			name="_data";
			break;
		case fitSet::fitToys:
			break;
	}
	name+="_";
	name+=_jptmin;
	name+="-";
	name+=_jptmax;
	name+="_";
	name+=_ymin;
	name+="-";
	name+=_ymax;
	if(_useEffs) {
		name+="_effWeight";
		if(_flavour==4) name+="_fit4";
		if(_flavour==5) name+="_fit5";
	}

	RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data"));
	const RooArgSet* vars = data->get();
	RooDataSet* ds(0);
	if(_useEffs) ds = new RooDataSet("ds", "ds", data, *vars, cutStr, weightName);
	else ds = new RooDataSet("ds", "ds", data, *vars, cutStr);
	//RooDataSet* ds = dynamic_cast<RooDataSet*>(ws->data("data")->reduce(cutStr));
	name+="_";
	name+=_ptBinsD[ptBinMin];
	name+="-";
	name+=_ptBinsD[ptBinMax];

	RooAbsPdf* data_pdf(0);

	TString toFloatStr, toFixStr, toZeroStr, toOneStr;
	switch(whichFit) {
		case fitType::fitMass:
			data_pdf = ws->pdf("sim_mass_pdf");
			typeName = "Mass";
			switch(whichSet) {
				case fitSet::fitJustCharm:
				case fitSet::fitJustBeauty:
				case fitSet::fitJustSignal:
					toFloatStr="signalYield*,dMean*,dBaseWidth*,dRatio*,dAlpha*,dFrac*";
					toFixStr="p0*,dN*";
					toZeroStr="bkgrndYield*";
					toOneStr="dWidthScale*";
					break;
				case fitSet::fitJustSideband:
					toFloatStr="bkgrndYield*,p0*";
					toFixStr="dMean*,dBaseWidth*,dWidthScale*,dRatio*,dAlpha*,dFrac*,dN*";
					toZeroStr="signalYield*";
					toOneStr="";
					break;
				case fitSet::fitData:
				case fitSet::fitToys:
					toFloatStr="signalYield*,bkgrndYield*,dMean*,dWidthScale*,p0*";
					toFixStr="dBaseWidth*,dRatio*,dAlpha*,dFrac*,dN*";
					toZeroStr="";
					toOneStr="";
					break;
			}
			break;
		case fitType::fitIP:
			data_pdf = ws->pdf("sim_ip_pdf");
			typeName = "IP";
			switch(whichSet) {
				case fitSet::fitJustCharm:
					toFloatStr="promptYield*,promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*";
					toFixStr="promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="displcYield*,bkgrndYield*,promptMeanShift*";
					toOneStr="";
					break;
				case fitSet::fitJustBeauty:
					toFloatStr="displcYield*,displcMean*,displcWidth*,displcAsym*";
					toFixStr="promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcRho*";
					toZeroStr="promptYield*,bkgrndYield*";
					toOneStr="";
					break;
				case fitSet::fitJustSignal:
					toFloatStr="promptYield*,displcYield*";
					toFixStr="promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="bkgrndYield*,promptMeanShift*";
					toOneStr="";
					break;
				case fitSet::fitJustSideband:
					toFloatStr="bkgrndYield*";
					toFixStr="promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="promptYield*,displcYield*";
					toOneStr="";
					break;
				case fitSet::fitData:
				case fitSet::fitToys:
					toFloatStr="promptYield*,displcYield*,bkgrndYield*,promptMeanShift*";
					toFixStr="promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="";
					toOneStr="";
					break;
			}
			break;
		case fitType::fit2D:
			data_pdf = ws->pdf("sim_data_pdf");
			typeName = "2D";
			switch(whichSet) {
				case fitSet::fitJustCharm:
					toFloatStr="promptYield*,dMean*,dBaseWidth*,promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*";
					toFixStr="dRatio*,dAlpha*,dFrac,dN*,*p0*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="displcYield*,bkgrndYield*,promptMeanShift*";
					toOneStr="dWidthScale*";
					break;
				case fitSet::fitJustBeauty:
					toFloatStr="displcYield*,dMean*,dBaseWidth*,displcMean*,displcWidth*,displcAsym*";
					toFixStr="dRatio*,dAlpha*,dFrac,dN*,p0*,promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcRho*";
					toZeroStr="promptYield*,bkgrndYield*";
					toOneStr="dWidthScale*";
					break;
				case fitSet::fitJustSignal:
					toFloatStr="promptYield*,displcYield*,dMean*,dBaseWidth*";
					toFixStr="dRatio*,dAlpha*,dFrac,dN*,p0*,promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="bkgrndYield*,promptMeanShift*";
					toOneStr="dWidthScale*";
					break;
				case fitSet::fitJustSideband:
					toFloatStr="bkgrndYield*,p0*";
					toFixStr="dMean*,dBaseWidth*,dWidthScale*,dRatio*,dAlpha*,dFrac*,dN*,promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="promptYield*,displcYield*";
					toOneStr="";
					break;
				case fitSet::fitData:
				case fitSet::fitToys:
					toFloatStr="promptYield*,displcYield*,bkgrndYield*,dMean*,dWidthScale*,p0*,promptMeanShift*";
					toFixStr="dBaseWidth*,dRatio*,dAlpha*,dFrac*,dN*,promptBaseMean*,promptWidth*,promptAsym*,promptRhoL*,promptRhoR*,displcMean*,displcWidth*,displcAsym*,displcRho*";
					toZeroStr="";
					toOneStr="";
					break;
			}
			break;
	}

	//unused PT bins
	for(uint ibin=0; ibin<_ptBinsD.size()-1; ++ibin) {
		if(ibin<ptBinMin || ibin>=ptBinMax) {
			toFixStr+=",*PTbin";
			toFixStr+=ibin;
			toZeroStr+=",*Yield*PTbin";
			toZeroStr+=ibin;
		}
	}

	//std::cout << toFloatStr << std::endl;//TODO
	//std::cout << toFixStr << std::endl;//TODO
	//std::cout << toZeroStr << std::endl;//TODO

	RooArgSet*  allPars   = data_pdf->getParameters(*ds);
	RooArgList* floatPars = static_cast<RooArgList*>(allPars->selectByName(toFloatStr));
	RooArgList* fixPars   = static_cast<RooArgList*>(allPars->selectByName(toFixStr));
	RooArgList* zeroPars  = static_cast<RooArgList*>(allPars->selectByName(toZeroStr));
	RooArgList* onePars   = static_cast<RooArgList*>(allPars->selectByName(toOneStr));
	RooArgList* yield4Pars= static_cast<RooArgList*>(allPars->selectByName("*promptYield*"));
	RooArgList* yield5Pars= static_cast<RooArgList*>(allPars->selectByName("*displcYield*"));

	for ( int i=0; i<floatPars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(floatPars->at(i))->setConstant(false);
	}
	for ( int i=0; i<fixPars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(fixPars->at(i))->setConstant(true);
	}
	for ( int i=0; i<onePars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(onePars->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(onePars->at(i))->setVal(1.0);
	}
	for ( int i=0; i<zeroPars->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(zeroPars->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(zeroPars->at(i))->setVal(0.0);
	}

	//check for unused or overlapping parameters
	RooArgList missingPars(*allPars);
	missingPars.remove(*floatPars);
	missingPars.remove(*fixPars);
	missingPars.remove(*zeroPars);
	missingPars.remove(*onePars);
	if(!missingPars.empty()) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were not specified" << std::endl;
		missingPars.printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*fixPars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be both floated and fixed" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*fixPars))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*zeroPars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be both floated and fixed to zero" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*zeroPars))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*onePars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be both floated and fixed to one" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*onePars))->printMultiline(std::cout,0xFF);
	}
	if(zeroPars->overlaps(*onePars)) {
		std::cout << "WARNING in SimDFitter::fit: the following parameters were specified to be fixed to both zero and one" << std::endl;
		static_cast<RooArgList*>(zeroPars->selectCommon(*onePars))->printMultiline(std::cout,0xFF);
	}

	//set sensible starting values for yields
	for(uint ibin=ptBinMin; ibin<ptBinMax; ++ibin) {
		//extract the yield in each bin from the dataset
		TString ptCutStr="ptCat==";
		ptCutStr+=ibin;
		double binYield = ds->sumEntries(ptCutStr);
		//std::cout << "entries" << ibin << " " << binYield << std::endl;

		if(binYield>=1) {
			TString yieldNameStr="*Yield*PTbin";
			yieldNameStr+=ibin;
			RooArgList* yields= static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(yieldNameStr)->selectByAttrib("Constant",kFALSE));
			for ( int i=0; i<yields->getSize(); ++i) {
				dynamic_cast<RooRealVar*>(yields->at(i))->setVal(binYield/yields->getSize());
				//std::cout << dynamic_cast<RooRealVar*>(yields->at(i))->getTitle() << " set to " << dynamic_cast<RooRealVar*>(yields->at(i))->getValV() << std::endl;
			}
		} else {
			std::cout << "Bin " << ibin << " has yield of " << binYield << ". Fixing to 0..." << std::endl;
			toFixStr+=",*PTbin";
			toFixStr+=ibin;
			toZeroStr+=",*Yield*PTbin";
			toZeroStr+=ibin;
		}
	}

	RooArgList* fixPars2   = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(toFixStr));
	RooArgList* zeroPars2  = static_cast<RooArgList*>(data_pdf->getParameters(*ds)->selectByName(toZeroStr));

	if(floatPars->overlaps(*fixPars2)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to fixed due to low bin populations" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*fixPars2))->printMultiline(std::cout,0xFF);
	}
	if(floatPars->overlaps(*zeroPars2)) {
		std::cout << "INFO in SimDFitter::fit: the following parameters were overridden from floated to zeroed due to low bin populations" << std::endl;
		static_cast<RooArgList*>(floatPars->selectCommon(*zeroPars2))->printMultiline(std::cout,0xFF);
	}

	//repeat in case more have been added
	for ( int i=0; i<fixPars2->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(fixPars2->at(i))->setConstant(true);
	}
	for ( int i=0; i<zeroPars2->getSize(); ++i) {
		dynamic_cast<RooRealVar*>(zeroPars2->at(i))->setConstant(true);
		dynamic_cast<RooRealVar*>(zeroPars2->at(i))->setVal(0.0);
	}

	//Print final list of floated parameters before fitting
	//std::cout << "Fit parameters:" << std::endl;
	//static_cast<RooArgList*>(allPars->selectByAttrib("Constant",kFALSE))->printMultiline(std::cout, 0xFF);
	
	double minCovQual(3);
	minCovQual=2;//allow not pos def//TODO
	gSystem->RedirectOutput(gSaveDir+"/log/fitFull"+typeName+name+".log","w");
	RooAbsReal* nll = data_pdf->createNLL(*ds, RooFit::Extended(), RooFit::NumCPU(4), RooFit::Range("FIT"));
	RooMinuit m(*nll);
	m.setPrintEvalErrors(-1);
	std::cout << "FIT STAGE 1" << std::endl;//TODO
	m.migrad();
	RooFitResult* r1 = m.save();
	std::cout << "STATUS=" << r1->status() << ", COVQUAL=" << r1->covQual() << std::endl;//TODO

	std::cout << "FIT STAGE 2" << std::endl;//TODO
	m.hesse();
	RooFitResult* r2 = m.save();
	std::cout << "STATUS=" << r2->status() << ", COVQUAL=" << r2->covQual() << std::endl;//TODO
	RooFitResult *r3(0), *r4(0), *r5(0), *r6(0);
	if(r2->status()!=0 || r2->covQual()<minCovQual) {
		std::cout << "Hesse failed STATUS=" << r2->status() << ", COVQUAL=" << r2->covQual() << std::endl;
		std::cout << "Retrying fit" << std::endl;
		std::cout << "FIT STAGE 3" << std::endl;//TODO
		m.migrad();
		r3 = m.save();
		std::cout << "STATUS=" << r3->status() << ", COVQUAL=" << r3->covQual() << std::endl;//TODO
		std::cout << "FIT STAGE 4" << std::endl;//TODO
		m.hesse();
		r4 = m.save();
		std::cout << "STATUS=" << r4->status() << ", COVQUAL=" << r4->covQual() << std::endl;//TODO
		if(r4->status()!=0 || r4->covQual()<minCovQual) {
			std::cout << "Hesse failed STATUS=" << r4->status() << ", COVQUAL=" << r4->covQual() << std::endl;
			//if(r4->status()!=0 || r4->covQual()<2) {//TODO allow covQual==2 for now
			status=false;
			//}
		}
	}
	//if we're using efficiency weights then re-run the fit using weights to calculate uncertainties
	if(_useEffs && !_skipSumW2Fits) {
		std::cout << "FIT STAGE 5" << std::endl;//TODO
		RooFitResult* r5 = data_pdf->fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(true), RooFit::PrintLevel(3));//TODO
		std::cout << "STATUS=" << r5->status() << ", COVQUAL=" << r5->covQual() << std::endl;//TODO
		status = (r5->status()==0);
		if(r5->status()!=0) {
			std::cout << "FIT STAGE 6" << std::endl;//TODO
			r6 = data_pdf->fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(true));
			std::cout << "STATUS=" << r6->status() << ", COVQUAL=" << r6->covQual() << std::endl;//TODO
			status = (r6->status()==0);
		}
	}
	gSystem->RedirectOutput(0);
	//static_cast<RooArgList*>(allPars->selectByAttrib("Constant",kFALSE))->printLatex(std::cout, 1, "NEYU", 2);
	delete nll;
	std::ofstream logfile(gSaveDir+"/log/fit"+typeName+name+".log",std::ofstream::out);
	r1->printMultiline(logfile,0);
	r2->printMultiline(logfile,0);
	if(r3) r3->printMultiline(logfile,0);
	if(r4) r4->printMultiline(logfile,0);
	if(r5) r5->printMultiline(logfile,0);
	if(r6) r6->printMultiline(logfile,0);
	logfile.close();

	if(true) { //status) {
		if(whichSet==fitSet::fitToys) {
			//do toy fits
			const RooCategory* cat = dynamic_cast<const RooCategory*>(&dynamic_cast<RooSimultaneous*>(data_pdf)->indexCat());
			RooMCStudy mcs(*data_pdf, RooArgSet(*DM,*DIP,*cat), RooFit::Extended(), RooFit::FitOptions(RooFit::PrintEvalErrors(-1)));
			mcs.generateAndFit(nToys, ds->sumEntries());
			TCanvas c1;
			TString varName;
			RooRealVar* var(0);
			for(uint ibin=ptBinMin; ibin<ptBinMax; ++ibin) {
				varName = "promptYield_PTbin";
				varName+=ibin;
				var = ws->var(varName);
				RooPlot* promptYieldPlot = mcs.plotPull(*var,RooFit::FitGauss(),RooFit::FrameBins(50),RooFit::FrameRange(-5,5));
				promptYieldPlot->SetTitle("");
				promptYieldPlot->SetXTitle("Pull(N_{prompt})");
				promptYieldPlot->Draw();
				c1.SaveAs(gSaveDir+"/fig/pull_"+varName+".pdf");

				varName = "displcYield_PTbin";
				varName+=ibin;
				var = ws->var(varName);
				RooPlot* displYieldPlot  = mcs.plotPull(*var,RooFit::FitGauss(),RooFit::FrameBins(50),RooFit::FrameRange(-5,5));
				displYieldPlot->SetTitle("");
				displYieldPlot->SetXTitle("Pull(N_{displ.})");
				displYieldPlot->Draw();
				c1.SaveAs(gSaveDir+"/fig/pull_"+varName+".pdf");
			}
		}
		//gSystem->RedirectOutput(0);
		if(whichSet!=fitSet::fitToys) {
//			double yieldS = signalYield.getVal();
//			double yieldB = bkgrndYield.getVal();
//			std::cout << yieldS << "\t" << yieldB << std::endl;
			
			std::vector<std::string> sig_pdfs;
			if(whichFit==fitType::fitMass && whichSet!=fitSet::fitJustSideband) {
				sig_pdfs.push_back( "signal*" );
			}
			if(whichFit!=fitType::fitMass && whichSet!=fitSet::fitJustSideband && whichSet!=fitSet::fitJustBeauty) {
				sig_pdfs.push_back( "prompt*" );
			}
			if(whichFit!=fitType::fitMass && whichSet!=fitSet::fitJustSideband && whichSet!=fitSet::fitJustCharm) {
				sig_pdfs.push_back( "displc*" );
			}
			std::vector<std::string> bkg_pdfs;
			if(whichSet!=fitSet::fitJustCharm && whichSet!=fitSet::fitJustBeauty && whichSet!=fitSet::fitJustSignal) {
				bkg_pdfs.push_back( "bkgrnd*" );
			}

			//for single-bin fits just plot that bin
			//int plotBins(-1);
			//if(ptBinMax==ptBinMin+1) {
			//	plotBins=ptBinMin;
			//}

			switch(whichFit) {
				case fitType::fitMass:
					//plotFit(*DM,     1784., 1944., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fitMass"+name, "#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]", plotBins);
					makeMassPlots(ds,name,ptBinMin,ptBinMax);
					break;
				case fitType::fitIP:
					//plotFit(*DIP,      -5.,   15., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fitIP"+name, "log #chi_{IP}^{2}", plotBins);
					makeIPPlots(ds,name,ptBinMin,ptBinMax);
					break;
				case fitType::fit2D:
					//plotFit(*DM,     1784., 1944., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fit2D_Mass"+name, "#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]", plotBins);
					//plotFit(*DIP,      -5.,   15., 80, ds, *data_pdf, sig_pdfs, bkg_pdfs, "fig/fit2D_IP"+name, "log #chi_{IP}^{2}", plotBins);
					make2DPlots(ds,name,ptBinMin,ptBinMax);
					printAllParams("dat/paramsAll"+typeName+name+".dat",ds,data_pdf);
					break;
			}
		}

		if(whichFit==fitType::fit2D) {
			_yield4=0.;
			_error4=0.;
			_yield5=0.;
			_error5=0.;
			for ( int i=0; i<yield4Pars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(yield4Pars->at(i));
				//std::cout << v->getTitle() << " " << v->getVal() << " " << v->getError() << std::endl;//TODO
				_yield4+= v->getVal();
				_error4+= v->getError()*v->getError();
			}
			_error4 = TMath::Sqrt(_error4);
			//std::cout << "..." << std::endl;//TODO
			for ( int i=0; i<yield5Pars->getSize(); ++i) {
				RooRealVar* v = dynamic_cast<RooRealVar*>(yield5Pars->at(i));
				//std::cout << v->getTitle() << " " << v->getVal() << " " << v->getError() << std::endl;//TODO
				_yield5+= v->getVal();
				_error5+= v->getError()*v->getError();
			}
			_error5 = TMath::Sqrt(_error5);
		}

		printParams("dat/params"+typeName+name+".dat",ds,data_pdf);
	}

	delete ds;

	return status;
}

void SimDFitter::makeMassPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
	RooRealVar* DM = ws->var("D0M");

	uint nBins1D(80);
	uint nBins(160);

	double mmin(1784.), mmax(1944.);

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);

	//pT-integrated histograms
	TH1* dataHist   = static_cast<TH1*>(ds->createHistogram("dataHist",*DM,RooFit::Binning(nBins,mmin,mmax)));
	TH1* modelHist  = new TH1D("modelHist", "",nBins,mmin,mmax);
	TH1* signalHist = new TH1D("signalHist","",nBins,mmin,mmax);
	TH1* bkgrndHist = new TH1D("bkgrndHist","",nBins,mmin,mmax);

	//axis labels
	for(auto hist : {dataHist,modelHist,signalHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
		hist->GetYaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	signalHist->SetLineColor(kRed);
	bkgrndHist->SetLineColor(kBlue);

	//store integrated and pT slice histograms
	std::vector<TH1*> dataHists {dataHist};
	std::vector<TH1*> modelHists {modelHist};
	std::vector<TH1*> signalHists {signalHist};
	std::vector<TH1*> bkgrndHists {bkgrndHist};

	//load all pT slices
	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
		RooAbsPdf* signal = ws->pdf(TString("signalMass_PTbin")+ibin);
		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrndMass_PTbin")+ibin);

		//if no parameters are split then only a single PDF will exist
		if(!signal) signal = ws->pdf(TString("signalMass"));
		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrndMass"));

		TH1* signalSlice = static_cast<TH1*>(signal->createHistogram(TString::Format("signalSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax)));
		TH1* bkgrndSlice = static_cast<TH1*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax)));
		TH1* dataSlice   = static_cast<TH1*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));

		TH1* modelSlice  = new TH1D(TString::Format("modelSlice%d",ibin), "",nBins,mmin,mmax);

		signalSlice->Scale(ws->var(TString::Format("signalYield_PTbin%d",ibin))->getVal());
		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());

		//make composites
		modelSlice->Add(signalSlice);
		modelSlice->Add(bkgrndSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		signalHist->Add(signalSlice);
		bkgrndHist->Add(bkgrndSlice);

		//setup style
		for(auto hist : {dataSlice,modelSlice,signalSlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
			hist->GetYaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		signalSlice->SetLineColor(kRed);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		signalHists.push_back(signalSlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH1* dataHist50 = dataHists[i]->Rebin( nBins/nBins1D, "dataHist50");

		//Mass plots
		TH1* dataHistMass   = dataHist50;
		TH1* signalHistMass = signalHists[i];
		TH1* bkgrndHistMass = bkgrndHists[i];
		TH1* modelHistMass  = modelHists[i];

		signalHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		bkgrndHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		modelHistMass ->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());

		dataHistMass->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f MeV/#it{c}^{2})", (mmax-mmin)/nBins1D));
		bkgrndHistMass->SetLineStyle(kDashed);
		signalHistMass->SetFillColor(kRed);
		signalHistMass->SetFillStyle(3350);

		dataHistMass->Draw("E1 X0 P");
		if(signalHistMass->Integral()>0) signalHistMass->Draw("HIST C SAME");
		if(bkgrndHistMass->Integral()>0) bkgrndHistMass->Draw("HIST C SAME");
		modelHistMass->Draw( "HIST C SAME");
		dataHistMass->Draw("E1 X0 P SAME");
		if(i==0) can.SaveAs(gSaveDir+"/fig/fitMass"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fitMass"+name+"PTbin"+(i-1)+".pdf");
	}
}

void SimDFitter::makeIPPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
	RooRealVar* DIP = ws->var("D0LogIPChi2");

	uint nBins1D(80);
	uint nBins(160);

	double ipmin(-5.), ipmax(15.);

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);

	//pT-integrated histograms
	TH1* dataHist   = static_cast<TH1*>(ds->createHistogram("dataHist",*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
	TH1* modelHist  = new TH1D("modelHist", "",nBins,ipmin,ipmax);
	TH1* signalHist = new TH1D("signalHist","",nBins,ipmin,ipmax);
	TH1* promptHist = new TH1D("promptHist","",nBins,ipmin,ipmax);
	TH1* displcHist = new TH1D("displcHist","",nBins,ipmin,ipmax);
	TH1* bkgrndHist = new TH1D("bkgrndHist","",nBins,ipmin,ipmax);

	//axis labels
	for(auto hist : {dataHist,modelHist,signalHist,promptHist,displcHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("log #chi_{IP}^{2}");
		hist->GetYaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	signalHist->SetLineColor(kRed);
	promptHist->SetLineColor(kGreen+2);
	displcHist->SetLineColor(kMagenta);
	bkgrndHist->SetLineColor(kBlue);

	//store integrated and pT slice histograms
	std::vector<TH1*> dataHists {dataHist};
	std::vector<TH1*> modelHists {modelHist};
	std::vector<TH1*> signalHists {signalHist};
	std::vector<TH1*> promptHists {promptHist};
	std::vector<TH1*> displcHists {displcHist};
	std::vector<TH1*> bkgrndHists {bkgrndHist};

	//load all pT slices
	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
		RooAbsPdf* prompt = ws->pdf(TString("promptIP_PTbin")+ibin);
		RooAbsPdf* displc = ws->pdf(TString("displcIP_PTbin")+ibin);
		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrndIP_PTbin")+ibin);

		//if no parameters are split then only a single PDF will exist
		if(!prompt) prompt = ws->pdf(TString("promptIP"));
		if(!displc) displc = ws->pdf(TString("displcIP"));
		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrndIP"));

		TH1* promptSlice = static_cast<TH2*>(prompt->createHistogram(TString::Format("promptSlice%d",ibin),*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
		TH1* displcSlice = static_cast<TH2*>(displc->createHistogram(TString::Format("displcSlice%d",ibin),*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
		TH1* bkgrndSlice = static_cast<TH2*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DIP,RooFit::Binning(nBins,ipmin,ipmax)));
		TH1* dataSlice   = static_cast<TH2*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DIP,RooFit::Binning(nBins,ipmin,ipmax),
					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));

		TH1* modelSlice  = new TH1D(TString::Format("modelSlice%d",ibin), "",nBins,ipmin,ipmax);
		TH1* signalSlice = new TH1D(TString::Format("signalSlice%d",ibin),"",nBins,ipmin,ipmax);

		promptSlice->Scale(ws->var(TString::Format("promptYield_PTbin%d",ibin))->getVal());
		displcSlice->Scale(ws->var(TString::Format("displcYield_PTbin%d",ibin))->getVal());
		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());

		//make composites
		modelSlice->Add(promptSlice);
		modelSlice->Add(displcSlice);
		modelSlice->Add(bkgrndSlice);
		signalSlice->Add(promptSlice);
		signalSlice->Add(displcSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		signalHist->Add(signalSlice);
		promptHist->Add(promptSlice);
		displcHist->Add(displcSlice);
		bkgrndHist->Add(bkgrndSlice);

		//setup style
		for(auto hist : {dataSlice,modelSlice,signalSlice,promptSlice,displcSlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("log #chi_{IP}^{2}");
			hist->GetYaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		signalSlice->SetLineColor(kRed);
		promptSlice->SetLineColor(kGreen+2);
		displcSlice->SetLineColor(kMagenta);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		signalHists.push_back(signalSlice);
		promptHists.push_back(promptSlice);
		displcHists.push_back(displcSlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH1* dataHist50 = dataHists[i]->Rebin( nBins/nBins1D,"dataHist50");

		//IP plots
		TH1* dataHistIP   = dataHist50;
		TH1* signalHistIP = signalHists[i];
		TH1* promptHistIP = promptHists[i];
		TH1* displcHistIP = displcHists[i];
		TH1* bkgrndHistIP = bkgrndHists[i];
		TH1* modelHistIP  = modelHists[i];

		signalHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		promptHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		displcHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		bkgrndHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		modelHistIP ->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());

		dataHistIP->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f)", (ipmax-ipmin)/nBins1D));
		bkgrndHistIP->SetLineStyle(kDashed);
		signalHistIP->SetLineStyle(kDotted);
		promptHistIP->SetFillColor(kGreen+2);
		displcHistIP->SetFillColor(kMagenta);
		promptHistIP->SetFillStyle(3345);
		displcHistIP->SetFillStyle(3354);

		dataHistIP->Draw("E1 X0 P");
		if(promptHistIP->Integral()>0 && displcHistIP->Integral()>0) signalHistIP->Draw("HIST C SAME");
		if(promptHistIP->Integral()>0) promptHistIP->Draw("HIST C SAME");
		if(displcHistIP->Integral()>0) displcHistIP->Draw("HIST C SAME");
		if(bkgrndHistIP->Integral()>0) bkgrndHistIP->Draw("HIST C SAME");
		modelHistIP->Draw( "HIST C SAME");
		dataHistIP->Draw("E1 X0 P SAME");
		if(i==0) can.SaveAs(gSaveDir+"/fig/fitIP"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fitIP"+name+"PTbin"+(i-1)+".pdf");
	}
}

void SimDFitter::make2DPlots(RooAbsData* ds, TString name, int minPTBin, int maxPTBin) {
	RooRealVar* DM = ws->var("D0M");
	RooRealVar* DIP = ws->var("D0LogIPChi2");

	uint nBins1D(50);
	uint nBins2D(30);
	uint nBins(150);

	double mmin(1834.), mmax(1894.), ipmin(-5.), ipmax(10.);

	//plotting style
	setLHCbStyle();
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
	Double_t greens[NRGBs] = { 1.00, 0.95, 0.50, 0.00, 0.00};
	Double_t blues[NRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);

	//pT-integrated histograms
	TH2* dataHist   = static_cast<TH2*>(ds->createHistogram("dataHist",*DM,RooFit::Binning(nBins,mmin,mmax),RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
	TH2* modelHist  = new TH2D("modelHist", "",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* signalHist = new TH2D("signalHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* promptHist = new TH2D("promptHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* displcHist = new TH2D("displcHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);
	TH2* bkgrndHist = new TH2D("bkgrndHist","",nBins,mmin,mmax,nBins,ipmin,ipmax);

	//axis labels
	for(auto hist : {dataHist,modelHist,signalHist,promptHist,displcHist,bkgrndHist}) {
		hist->SetTitle("");
		hist->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
		hist->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
		hist->GetZaxis()->SetTitle("");
	}

	dataHist->SetLineColor(kBlack);
	dataHist->SetMarkerStyle(20);
	modelHist->SetLineColor(kBlack);
	signalHist->SetLineColor(kRed);
	promptHist->SetLineColor(kGreen+2);
	displcHist->SetLineColor(kMagenta);
	bkgrndHist->SetLineColor(kBlue);
	promptHist->SetContour(10);
	displcHist->SetContour(10);
	bkgrndHist->SetContour(5);

	//store integrated and pT slice histograms
	std::vector<TH2*> dataHists {dataHist};
	std::vector<TH2*> modelHists {modelHist};
	std::vector<TH2*> signalHists {signalHist};
	std::vector<TH2*> promptHists {promptHist};
	std::vector<TH2*> displcHists {displcHist};
	std::vector<TH2*> bkgrndHists {bkgrndHist};

	//load all pT slices
	for(int ibin=minPTBin; ibin<maxPTBin; ++ibin) {
		RooAbsPdf* prompt = ws->pdf(TString("prompt_PTbin")+ibin);
		RooAbsPdf* displc = ws->pdf(TString("displc_PTbin")+ibin);
		RooAbsPdf* bkgrnd = ws->pdf(TString("bkgrnd_PTbin")+ibin);

		//if no parameters are split then only a single PDF will exist
		if(!prompt) prompt = ws->pdf(TString("prompt"));
		if(!displc) displc = ws->pdf(TString("displc"));
		if(!bkgrnd) bkgrnd = ws->pdf(TString("bkgrnd"));

		TH2* promptSlice = static_cast<TH2*>(prompt->createHistogram(TString::Format("promptSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
		TH2* displcSlice = static_cast<TH2*>(displc->createHistogram(TString::Format("displcSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
		TH2* bkgrndSlice = static_cast<TH2*>(bkgrnd->createHistogram(TString::Format("bkgrndSlice%d",ibin),*DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax))));
		TH2* dataSlice   = static_cast<TH2*>(ds->createHistogram(    TString::Format("dataSlice%d",ibin),  *DM,RooFit::Binning(nBins,mmin,mmax),
					                                                             RooFit::YVar(*DIP,RooFit::Binning(nBins,ipmin,ipmax)),
					                                                             RooFit::Cut(TString::Format("ptCat==%d",ibin))));

		TH2* modelSlice  = new TH2D(TString::Format("modelSlice%d",ibin), "",nBins,mmin,mmax,nBins,ipmin,ipmax);
		TH2* signalSlice = new TH2D(TString::Format("signalSlice%d",ibin),"",nBins,mmin,mmax,nBins,ipmin,ipmax);

		promptSlice->Scale(ws->var(TString::Format("promptYield_PTbin%d",ibin))->getVal());
		displcSlice->Scale(ws->var(TString::Format("displcYield_PTbin%d",ibin))->getVal());
		bkgrndSlice->Scale(ws->var(TString::Format("bkgrndYield_PTbin%d",ibin))->getVal());

		//make composites
		modelSlice->Add(promptSlice);
		modelSlice->Add(displcSlice);
		modelSlice->Add(bkgrndSlice);
		signalSlice->Add(promptSlice);
		signalSlice->Add(displcSlice);

		//add to integrated histograms
		modelHist->Add(modelSlice);
		signalHist->Add(signalSlice);
		promptHist->Add(promptSlice);
		displcHist->Add(displcSlice);
		bkgrndHist->Add(bkgrndSlice);

		//setup style
		for(auto hist : {dataSlice,modelSlice,signalSlice,promptSlice,displcSlice,bkgrndSlice}) {
			hist->SetTitle("");
			hist->GetXaxis()->SetTitle("#it{m}(#it{K}^{#minus}#pi^{#plus}) [MeV/#it{c}^{2}]");
			hist->GetYaxis()->SetTitle("log #chi_{IP}^{2}");
			hist->GetZaxis()->SetTitle("");
		}

		dataSlice->SetLineColor(kBlack);
		dataSlice->SetMarkerStyle(20);
		modelSlice->SetLineColor(kBlack);
		signalSlice->SetLineColor(kRed);
		promptSlice->SetLineColor(kGreen+2);
		displcSlice->SetLineColor(kMagenta);
		bkgrndSlice->SetLineColor(kBlue);

		dataHists.push_back(dataSlice);
		modelHists.push_back(modelSlice);
		signalHists.push_back(signalSlice);
		promptHists.push_back(promptSlice);
		displcHists.push_back(displcSlice);
		bkgrndHists.push_back(bkgrndSlice);
	}

	//make pull histograms
	TH2* dataHist10 = dataHist->Rebin2D(nBins/10,nBins/10,"dataHist10");
	TH2* modelHist10 = modelHist->Rebin2D(nBins/10,nBins/10,"modelHist10");

	//create pull histogram
	modelHist10->Scale(dataHist10->Integral()/modelHist10->Integral());
	for(int ibin=1; ibin<=10; ++ibin) {
		for(int jbin=1; jbin<=10; ++jbin) {
			double D = dataHist10->GetBinContent(ibin,jbin);
			double F = modelHist10->GetBinContent(ibin,jbin);
			double E = dataHist10->GetBinError(ibin,jbin);
			modelHist10->SetBinContent(ibin,jbin,(D-F)/E);
		}
	}

	TCanvas can;

	for(uint i=0; i<dataHists.size(); ++i) {
		TH2* dataHist50 = dataHists[i]->Rebin2D( nBins/nBins1D, nBins/nBins1D,"dataHist50");

		//Mass plots
		TH1* dataHistMass   = dataHist50->ProjectionX("dataHistMass");
		TH1* signalHistMass = signalHists[i]->ProjectionX("signalHistMass");
		TH1* promptHistMass = promptHists[i]->ProjectionX("promptHistMass");
		TH1* displcHistMass = displcHists[i]->ProjectionX("displcHistMass");
		TH1* bkgrndHistMass = bkgrndHists[i]->ProjectionX("bkgrndHistMass");
		TH1* modelHistMass  = modelHists[i]->ProjectionX("modelHistMass");

		signalHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		promptHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		displcHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		bkgrndHistMass->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());
		modelHistMass ->Scale(nBins/nBins1D*dataHistMass->Integral()/modelHistMass->Integral());

		dataHistMass->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f MeV/#it{c}^{2})", (mmax-mmin)/nBins1D));
		bkgrndHistMass->SetLineStyle(kDashed);
		signalHistMass->SetLineStyle(kDotted);
		promptHistMass->SetFillColor(kGreen+2);
		displcHistMass->SetFillColor(kMagenta);
		promptHistMass->SetFillStyle(3345);
		displcHistMass->SetFillStyle(3354);

		dataHistMass->Draw("E1 X0 P");
		signalHistMass->Draw("HIST C SAME");
		promptHistMass->Draw("HIST C SAME");
		displcHistMass->Draw("HIST C SAME");
		bkgrndHistMass->Draw("HIST C SAME");
		modelHistMass->Draw( "HIST C SAME");
		dataHistMass->Draw("E1 X0 P SAME");
		if(i==0) can.SaveAs(gSaveDir+"/fig/fit2D_Mass"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fit2D_Mass"+name+"PTbin"+(i-1)+".pdf");

		//IP plots
		TH1* dataHistIP   = dataHist50->ProjectionY("dataHistIP");
		TH1* signalHistIP = signalHists[i]->ProjectionY("signalHistIP");
		TH1* promptHistIP = promptHists[i]->ProjectionY("promptHistIP");
		TH1* displcHistIP = displcHists[i]->ProjectionY("displcHistIP");
		TH1* bkgrndHistIP = bkgrndHists[i]->ProjectionY("bkgrndHistIP");
		TH1* modelHistIP  = modelHists[i]->ProjectionY("modelHistIP");

		signalHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		promptHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		displcHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		bkgrndHistIP->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());
		modelHistIP ->Scale(nBins/nBins1D*dataHistIP->Integral()/modelHistIP->Integral());

		dataHistIP->GetYaxis()->SetTitle(TString::Format("Candidates / (%.2f)", (ipmax-ipmin)/nBins1D));
		bkgrndHistIP->SetLineStyle(kDashed);
		signalHistIP->SetLineStyle(kDotted);
		promptHistIP->SetFillColor(kGreen+2);
		displcHistIP->SetFillColor(kMagenta);
		promptHistIP->SetFillStyle(3345);
		displcHistIP->SetFillStyle(3354);

		dataHistIP->Draw("E1 X0 P");
		signalHistIP->Draw("HIST C SAME");
		promptHistIP->Draw("HIST C SAME");
		displcHistIP->Draw("HIST C SAME");
		bkgrndHistIP->Draw("HIST C SAME");
		modelHistIP->Draw( "HIST C SAME");
		dataHistIP->Draw("E1 X0 P SAME");
		if(i==0) can.SaveAs(gSaveDir+"/fig/fit2D_IP"+name+".pdf");
		else can.SaveAs(gSaveDir+"/fig/fit2D_IP"+name+"PTbin"+(i-1)+".pdf");
	}
	
	//2D plots
	TH2* dataHist30 = dataHist->Rebin2D( nBins/nBins2D, nBins/nBins2D,"dataHist30");
	can.SetTickx();
	can.SetTicky();
	dataHist30->Draw("col");
	//signalHist->Draw("cont2 same");
	promptHist->Draw("cont2 same");
	displcHist->Draw("cont2 same");
	bkgrndHist->Draw("cont2 same");
	can.SaveAs(gSaveDir+"/fig/data"+name+".pdf");
	promptHist->Draw("col");
	can.SaveAs(gSaveDir+"/fig/promptPDF"+name+".pdf");
	displcHist->Draw("col");
	can.SaveAs(gSaveDir+"/fig/displcPDF"+name+".pdf");
	bkgrndHist->Draw("col");
	can.SaveAs(gSaveDir+"/fig/bkgrndPDF"+name+".pdf");

	can.SetRightMargin(0.10);
	const Int_t NRGBs2 = 4;
	const Int_t NCont2 = 255;
	Double_t stops2[NRGBs2]  = { 0.00, 0.35, 0.65, 1.00};
	Double_t reds2[NRGBs2]   = { 0.00, 1.00, 1.00, 1.00};
	Double_t greens2[NRGBs2] = { 0.00, 1.00, 1.00, 0.00};
	Double_t blues2[NRGBs2]  = { 1.00, 1.00, 1.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);
	gStyle->SetNumberContours(NCont2);
	modelHist10->SetMinimum(-3.);
	modelHist10->SetMaximum( 3.);
	modelHist10->Draw("colz");
	can.SaveAs(gSaveDir+"/fig/pull"+name+".pdf");
}

bool SimDFitter::doFits() {
	if(!_initialised) init();
	ws->writeToFile("testWS.root");

	bool status(true);
	int nRetries(0);
	int nAttemptsMC(1);
	int nAttemptsData(2);
	if(!_mcFitsDone) {
		do {
			std::cout << "Fitting mass: simulated signal" << std::endl;
			status=fit(fitType::fitMass, fitSet::fitJustSignal,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsMC);
		//do {
		//	std::cout << "Fitting mass: simulated background" << std::endl;
		//	status=fit(fitType::fitMass, fitSet::fitJustSideband,0,_ptBinsD.size()-1);
		//	if(status) break;
		//	std::cout << "Failed: retrying..." << std::endl;
		//	++nRetries;
		//} while(!status && nRetries<nAttemptsMC);

		if(!status) {
			std::cout << "mass fit failed, skipping IP fit" << std::endl;
			return status;
		}

		do {
			std::cout << "Fitting IP: simulated prompt" << std::endl;
			//for(uint ibin=0; ibin<_ptBinsD.size()-1; ++ibin) {
			//	status=fit(fitType::fitIP, fitSet::fitJustCharm,ibin,ibin+1);
			//	std::cout << ibin << " " << status << std::endl;//TODO
			//}
			status=fit(fitType::fitIP, fitSet::fitJustCharm,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsMC);
		if(!status) {
			std::cout << "IP fit failed, skipping further fits" << std::endl;
			return status;
		}
		do {
			std::cout << "Fitting IP: simulated displaced" << std::endl;
			status=fit(fitType::fitIP, fitSet::fitJustBeauty,0,_ptBinsD.size()-1);
			if(status) break;
			std::cout << "Failed: retrying..." << std::endl;
			++nRetries;
		} while(!status && nRetries<nAttemptsMC);
		if(!status) {
			std::cout << "IP fit failed, skipping further fits" << std::endl;
			return status;
		}
		//do {
		//	std::cout << "Fitting IP: simulated background" << std::endl;
		//	status=fit(fitType::fitIP, fitSet::fitJustSideband,0,_ptBinsD.size()-1);
		//	if(status) break;
		//	std::cout << "Failed: retrying..." << std::endl;
		//	++nRetries;
		//} while(!status && nRetries<nAttemptsMC);
		if(!status) {
			std::cout << "IP fit failed, skipping 2D fit" << std::endl;
			return status;
		}
		_mcFitsDone=true;
	}
	nRetries=0;//reset allowed retries for fit to data
	do {
		std::cout << "Fitting 2D: data" << std::endl;
		status=fit(fitType::fit2D, fitSet::fitData,0,_ptBinsD.size()-1);
		if(status) break;
		std::cout << "Failed: retrying..." << nRetries << std::endl;
		++nRetries;
	} while(!status && nRetries<nAttemptsData);
	return status;//TODO
	status&=fit(fitType::fit2D, fitSet::fitToys,0,_ptBinsD.size()-1,100);
}

double SimDFitter::getAveWeight(double flavour) {
	TFile* f = TFile::Open(dFileName());
	if(!f) return 1.;

	TTree* t7 = static_cast<TTree*>(f->Get("T7"));
	if(!t7) return 1.;

	double mmin(1834.), mmax(1894.);
	double ipmin(-5.), ipmax(10.);

	TString cutStr = TString::Format("(JetPT>=%f && JetPT<%f && JetEta>=%f && JetEta<%f)",_jptmin,_jptmax,_ymin,_ymax);
	TString weightStr;

	if(flavour==4) {
		weightStr="(weight4)*";
		ipmax= 3;
	} else if(flavour==5) {
		weightStr="(weight5)*";
		ipmin= 3;
	}
	TH2D sum("sum","",1,mmin,mmax,1,ipmin,ipmax);
	TH2D count("count","",1,mmin,mmax,1,ipmin,ipmax);
	t7->Draw("D0LogIPChi2:D0M>>sum",weightStr+cutStr);
	t7->Draw("D0LogIPChi2:D0M>>count",cutStr);

	return sum.GetBinContent(1,1)/count.GetBinContent(1,1);
}

//int main() {
//	RooMsgService::instance().setSilentMode(true);
//	RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
//	gErrorIgnoreLevel = kError;
//	gSaveDir="simFit/";
//	SimDFitter s;
//	s.doFits();
//}
