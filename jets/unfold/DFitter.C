#include "DFitter.h"

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
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMsgService.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include "RooPromptShape.h"

#include <boost/progress.hpp>

#include "cloneHists.h"
#include "outputFunctions.h"
#include "DatasetManager.h"

void DFitter::setInputs(TString data, TString eff, TString acc, bool simple, bool isMC, bool useRhoZCorr) {
	_dataFile  = data;
	_effFile   = eff;
	_accFile   = acc;
	_useSimpleEffs = simple;
	_dataIsMC = isMC;
	_useRhoZCorr = useRhoZCorr;

	_inputsSet=true;
}

void DFitter::setDPtRange(double ptmin, double ptmax) {
	_dptmin   = ptmin;
	_dptmax   = ptmax;
}

bool DFitter::addEffsFull() {
	std::cout << "INFO : adding efficiencies to file " << _name << std::endl;
	//TFile* f0 = TFile::Open(_dataFile);

	TFile* f2 = TFile::Open(_effFile);
	TFile* f3 = TFile::Open(_accFile);

	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset(_dataFile);
	//TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	//if(!t0) return false;

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
	if(_dataIsMC) dm->setBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = dm->getEntries();

	TFile* fout = TFile::Open(dFileName(),"RECREATE");
	TTree* tout = new TTree("T","");

	double D0M(0.), D0PT(0.), D0Eta(0.), D0LogIPChi2(0.), weight4(0.), weight5(0.);
	double D0RhoSq(0.), D0Z(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetEta",       &JetEta);
	tout->Branch("D0M",          &D0M);
	tout->Branch("D0PT",         &D0PT);
	tout->Branch("D0Eta",        &D0Eta);
	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
	tout->Branch("weight4",      &weight4);
	tout->Branch("weight5",      &weight5);

	if(_dataIsMC) tout->Branch("JetTruePT", &JetTruePT);

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

			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;

			//check PID cuts
			if(!_dataIsMC && vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			//if(!dataIsMC && vD0PIPNNPI->at(s)<0.1 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //pion PID turned off

			D0M         = vD0M->at(s);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
			D0PT        = D0P.Pt();
			D0Eta       = D0P.Eta();
			D0RhoSq     = vD0X->at(s)*vD0X->at(s) + vD0Y->at(s)*vD0Y->at(s);
			D0Z         = vD0Z->at(s);

			//get efficiency corrections
			double effacc(0.), effrec(0.), effcor(1.), effsel4(0.), effsel5(0.), effpid(0.);

			if(D0PT>=100000.) D0PT=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(D0PT/1000.,D0Eta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(D0PT/1000.,D0Eta));
			if(_useRhoZCorr) effcor = hcor ->GetBinContent(hcor ->FindBin(D0RhoSq   ,D0Z  ));
			effsel4=      hsel4->GetBinContent(hsel4->FindBin(D0PT      ,D0Eta));
			effsel5=      hsel5->GetBinContent(hsel5->FindBin(D0PT      ,D0Eta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(D0PT      ,D0Eta));

			double eff4=effacc*effrec*effcor*effsel4*effpid;
			double eff5=effacc*effrec*effcor*effsel5*effpid;

			if(eff4<0.01 || eff5<0.01) {
//TODO//				std::cout << D0PT << "\t" << D0Eta << "\t" << effacc << "\t" << effrec << "\t" << effcor << "\t" << effsel4 << "\t" << effsel5 << "\t" << effpid << std::endl;
//TODO//				if(effcor==0) std::cout << D0RhoSq << "\t" << D0Z << std::endl;//TODO
				continue;
			}

			weight4 = 1./eff4;
			weight5 = 1./eff5;

			if(_dataIsMC) {
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
			break;//only keep one D0 candidate per entry
		}
	}

	dm->reset();

	tout->AutoSave();
	fout->Close();

	_effsSet = true;

	return true;
}

bool DFitter::addEffsSimple() {
	std::cout << "INFO : adding efficiencies to file " << _name << std::endl;
	//TFile* f0 = TFile::Open(_dataFile);

	TFile* f2 = TFile::Open(_effFile);

	DatasetManager* dm = DatasetManager::getInstance();
	dm->setDataset("data");
	//TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	//if(!t0) return false;

	TH2D* heff4 = dynamic_cast<TH2D*>(f2->Get("effD04"));
	TH2D* heff5 = dynamic_cast<TH2D*>(f2->Get("effD05"));

	if(!heff4 || !heff5) return false;

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

	dm->setBranchAddress("JetPT",         &JetPT);
	dm->setBranchAddress("JetEta",        &JetEta);
	if(_dataIsMC) dm->setBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = dm->getEntries();

	TFile* fout = TFile::Open(dFileName(),"RECREATE");
	TTree* tout = new TTree("T","");

	double D0M(0.), D0PT(0.), D0Eta(0.), D0LogIPChi2(0.), weight4(0.), weight5(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetEta",       &JetEta);
	tout->Branch("D0M",          &D0M);
	tout->Branch("D0PT",         &D0PT);
	tout->Branch("D0Eta",        &D0Eta);
	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
	tout->Branch("weight4",      &weight4);
	tout->Branch("weight5",      &weight5);

	if(_dataIsMC) tout->Branch("JetTruePT", &JetTruePT);

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

			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>_dptmin)) continue;
			if(!(_dptmax==-1 || D0P.Pt() < _dptmax) ) continue;

			//check PID cuts
			if(vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			//if(vD0PIPNNPI->at(s)<0.1) continue; //pion PID turned off

			D0M         = vD0M->at(s);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
			D0PT        = D0P.Pt();
			D0Eta       = D0P.Eta();

			//get efficiency corrections
			double eff4(0.), eff5(0.);

			if(D0PT>=100000.) D0PT=99999.;
			eff4   =      heff4->GetBinContent(heff4->FindBin(D0PT      ,D0Eta));
			eff5   =      heff5->GetBinContent(heff5->FindBin(D0PT      ,D0Eta));

			if(eff4<0.01 || eff5<0.01) {
//TODO//				std::cout << D0PT << "\t" << D0Eta << "\t" << eff4 << "\t" << eff5 << std::endl;
				continue;
			}

			weight4 = 1./eff4;
			weight5 = 1./eff5;

			//TODO//breaks eff if comparing against unweighted MC//if(_dataIsMC) {//scale weights for roughly continuous jet true pT
			//TODO//breaks eff if comparing against unweighted MC//	if(JetTruePT>50000.) {
			//TODO//breaks eff if comparing against unweighted MC//		weight4*=0.007;
			//TODO//breaks eff if comparing against unweighted MC//		weight5*=0.007;
			//TODO//breaks eff if comparing against unweighted MC//	} else if(JetTruePT>20000.) {
			//TODO//breaks eff if comparing against unweighted MC//		weight4*=0.10;
			//TODO//breaks eff if comparing against unweighted MC//		weight5*=0.10;
			//TODO//breaks eff if comparing against unweighted MC//	} else if(JetTruePT>15000.) {
			//TODO//breaks eff if comparing against unweighted MC//		weight4*=0.25;
			//TODO//breaks eff if comparing against unweighted MC//		weight5*=0.25;
			//TODO//breaks eff if comparing against unweighted MC//	}
			//TODO//breaks eff if comparing against unweighted MC//}

			tout->Fill();
			break;//only keep one D0 candidate per entry
		}
	}

	dm->reset();

	tout->AutoSave();
	fout->Close();

	_effsSet = true;

	return true;
}

bool DFitter::fitD_1D1D(double& yield, double& error) {
	std::cout << "INFO : fitting D0 - 1D mass and IP fits to file " << _name << ", pT in (" << _jptmin << "," << _jptmax << ")" << std::endl;
	bool fix(true);

	double valPromptMean(0.), valPromptWidth(0.), valPromptAsym(0.), valPromptRhoL(0.), valPromptRhoR(0.),  valDisplacedMean(5.), valDisplacedWidth(1.5);
	//D0 settings
	valPromptMean = 0.87;
	valPromptWidth = 1.09;
	valPromptAsym = -0.29;
	valPromptRhoL = 1.30;
	valPromptRhoR = 1.69;

	TString ptStr;
	ptStr+=_jptmin; ptStr+="-"; ptStr+=_jptmax;

	if(_dptmin!=5000. || _dptmax!=-1) {
		ptStr+="_D0";
		ptStr+=_dptmin;
		ptStr+="-";
		ptStr+=_dptmax;
	}
	if(!_useEffs) {
		ptStr+="_noEffCorr";
	}

	TFile* f0 = TFile::Open(dFileName());

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TString weightName;

	if(_flavour==4) {
		if(_useEffs) weightName = "weight4";
	} else if(_flavour==5) {
		if(_useEffs) weightName = "weight5";
	} else {
		return false;
	}

	TString dType("D0");
	double massVal(1865.);
	double massScaleFact(1.);
	double ratioVal(3.), fracVal(0.9), alphaVal(3.), nVal(1.);

	//D0 settings below
	massVal=1865.;
	ratioVal=2.85;
	fracVal=0.89;
	alphaVal=3.0;
	nVal=1.0;

	RooRealVar DM(  dType+"M",  dType+"M",  massVal-massScaleFact*80.,  massVal+massScaleFact*80.,  "MeV/#it{c}^{2}"); 
	DM.setBins(30);
	RooRealVar DLOGIPCHI2(  dType+"LogIPChi2",  dType+"LogIPChi2",  -5.,  15.,  ""); 
	DLOGIPCHI2.setBins(20);
	RooRealVar DPT(  dType+"PT",  dType+"PT",  0.,  1000000.); 
	RooRealVar JetPT( "JetPT",  "JetPT",  _jptmin,  _jptmax,  "MeV/#it{c}"); 
	RooRealVar DWeight(  weightName,  weightName,  0.,   1000.); 

	RooArgSet obs;
	obs.add(JetPT);
	obs.add(DM);
	obs.add(DPT);
	obs.add(DLOGIPCHI2);
	if(_useEffs) obs.add(DWeight);

	RooDataSet* ds(0);
	if(_useEffs) ds = new RooDataSet("ds","ds", obs, RooFit::WeightVar(weightName), RooFit::Import(*t0));
	else ds = new RooDataSet("ds","ds", obs, RooFit::Import(*t0));

	TString peakCut = "TMath::Abs("; peakCut+=dType; peakCut+="M-"; peakCut+=massVal; peakCut+=")<"; peakCut+=massScaleFact*20.;
	TString sideCut = "TMath::Abs("; sideCut+=dType; sideCut+="M-"; sideCut+=massVal; sideCut+=")<"; sideCut+=massScaleFact*80.; sideCut+=" && ";
	sideCut+= "TMath::Abs("; sideCut+=dType; sideCut+="M-"; sideCut+=massVal; sideCut+=")>"; sideCut+=massScaleFact*40.;

	RooDataSet* dsPeak = dynamic_cast<RooDataSet*>(ds->reduce(peakCut));
	RooDataSet* dsSB   = dynamic_cast<RooDataSet*>(ds->reduce(sideCut));

	//fit model parameters
	RooRealVar dMass("dMass","mass mean",massVal,massVal-massScaleFact*20.,massVal+massScaleFact*20.);
	RooRealVar dWidth("dWidth","mass width",massScaleFact*20.,massScaleFact*2.,massScaleFact*100.);

	RooRealVar dRatio("dRatio","ratio of widths", ratioVal, 1.0, 5);
	RooFormulaVar dWidthG("dWidthG","","@0*@1",RooArgList(dWidth,dRatio));
	RooRealVar dAlpha("dAlpha", "dAlpha", alphaVal, 0.5, 5.0);
	RooRealVar dN(    "dN",     "dN",     nVal);//2.0, 0.0, 5.0);
	RooRealVar dFrac("dFrac","fraction in CB", fracVal, 0., 1.);

	dRatio.setConstant();
	dFrac.setConstant();
	dAlpha.setConstant();
	dN.setConstant();

	RooGaussian sigMass_gauss("sigMass_gauss","",DM,dMass,dWidthG);
	RooCBShape  sigMass_cb("sigMass_cb", "", DM, dMass, dWidth, dAlpha, dN); 
	RooAddPdf sigMass( "sigMass", "", RooArgList(sigMass_cb,sigMass_gauss), RooArgList(dFrac) );

	RooRealVar p0("p0","p0",0., -0.1, 0.1);
	RooPolynomial bkgMass("bkgMass","",DM, RooArgList(p0));

	double maxyield = 1e4;
	double fracsig  = 0.9;
	TString cutStr="TMath::Abs(";
	cutStr+=dType;
	cutStr+="M-";
	cutStr+=massVal;
	cutStr+=")<";
	cutStr+=massScaleFact*80.;
	cutStr+=" && ";
	cutStr+=dType;
	cutStr+="LogIPChi2>-5. && ";
	cutStr+=dType;
	cutStr+="LogIPChi2<15.";
	TString cutStrPeak="TMath::Abs(";
	cutStrPeak+=dType;
	cutStrPeak+="M-";
	cutStrPeak+=massVal;
	cutStrPeak+=")<";
	cutStrPeak+=massScaleFact*30.;
	cutStrPeak+=" && ";
	cutStrPeak+=dType;
	cutStrPeak+="LogIPChi2>-5. && ";
	cutStrPeak+=dType;
	cutStrPeak+="LogIPChi2<15.";
	TString cutStrSB="TMath::Abs(";
	cutStrSB+=dType;
	cutStrSB+="M-";
	cutStrSB+=massVal;
	cutStrSB+=")<";
	cutStrSB+=massScaleFact*80.;
	cutStrSB+=" && ";
	cutStrSB+="TMath::Abs(";
	cutStrSB+=dType;
	cutStrSB+="M-";
	cutStrSB+=massVal;
	cutStrSB+=")>";
	cutStrSB+=massScaleFact*50.;
	cutStrSB+=" && ";
	cutStrSB+=dType;
	cutStrSB+="LogIPChi2>-5. && ";
	cutStrSB+=dType;
	cutStrSB+="LogIPChi2<15.";

	maxyield = t0->GetEntries(cutStr);
	fracsig = (t0->GetEntries(cutStrPeak)-t0->GetEntries(cutStrSB))/maxyield;
	std::cout << maxyield << "\t" << fracsig << std::endl;

	maxyield = ds->sumEntries();
	fracsig = (ds->sumEntries(cutStrPeak)-ds->sumEntries(cutStrSB))/maxyield;
	std::cout << maxyield << "\t" << fracsig << std::endl;

	// -- yields
	RooRealVar sigYield(  "sigYield",  "yield D",        fracsig *maxyield,    0.0,     1.1*maxyield);
	RooRealVar bkgYield(  "bkgYield",  "yield comb", (1.-fracsig)*maxyield,    0.0,     1.1*maxyield);

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(sigMass,bkgMass), RooArgList(sigYield,bkgYield) );

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/"+dType+"_"+_name+"_"+ptStr+typeName(_type)+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));
	gSystem->RedirectOutput(0);

	RooRealVar  bkgYield_mean("bkgYield_mean","bkgYield_mean", bkgYield.getVal());
	RooRealVar  bkgYield_sigma("bkgYield_sigma","bkgYield_sigma", bkgYield.getError());
	RooGaussian bkgYield_constraint("bkgYield_constraint","bkgYield_constraint", bkgYield, bkgYield_mean, bkgYield_sigma);

	DM.setRange("signal",massVal-massScaleFact*20.,massVal+massScaleFact*20.);
	DM.setRange("sideLo",massVal-massScaleFact*80.,massVal-massScaleFact*60.);
	DM.setRange("sideHi",massVal+massScaleFact*60.,massVal+massScaleFact*80.);
	DM.setRange("full",  massVal-massScaleFact*80.,massVal+massScaleFact*80.);

	DLOGIPCHI2.setRange("sigIPCHI2",-5.,3.);
	DLOGIPCHI2.setRange("fullIPCHI2",-5.,15.);

	double fsig_1 = sigMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("signal"))->getVal();
	double fsig_2 = sigMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("sideLo"))->getVal();
	double fsig_3 = sigMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("sideHi"))->getVal();
	double fsig_0 = sigMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("full"))->getVal();

	double fbkg_1 = bkgMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("signal"))->getVal();
	double fbkg_2 = bkgMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("sideLo"))->getVal();
	double fbkg_3 = bkgMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("sideHi"))->getVal();
	double fbkg_0 = bkgMass.createIntegral(RooArgSet(DM),RooFit::NormSet(DM),RooFit::Range("full"))->getVal();

	std::cout << "yields in mass windows" << std::endl;
	std::cout << sigYield.getVal()*fsig_1/fsig_0 << "\t" << bkgYield.getVal()*fbkg_1/fbkg_0 << std::endl;
	std::cout << sigYield.getVal()*(fsig_2+fsig_3)/fsig_0 << "\t" << bkgYield.getVal()*(fbkg_2+fbkg_3)/fbkg_0 << std::endl;

	RooRealVar promptMean("promptMean","mean prompt",valPromptMean,-1.,3.);
	RooRealVar promptWidth("promptWidth","width prompt",valPromptWidth,0.5,5.);
	RooRealVar promptAsym("promptAsym","asym. prompt",valPromptAsym,-1.,1.);
	RooRealVar promptRhoL("promptRhoL","exp L prompt",valPromptRhoL,0.5,3.);
	RooRealVar promptRhoR("promptRhoR","exp R prompt",valPromptRhoR,0.5,3.);
	RooPromptShape promptLOGIPCHI2("promptLOGIPCHI2","",DLOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	RooRealVar displacedMean("displacedMean","mean displ.",valDisplacedMean,2.,12.);
	RooRealVar displacedWidth("displacedWidth","width displ.",valDisplacedWidth,0.5,4.);
	RooGaussian displacedLOGIPCHI2("displacedLOGIPCHI2","",DLOGIPCHI2,displacedMean,displacedWidth);

	RooKeysPdf bkgLOGIPCHI2("bkgLOGIPCHI2","",DLOGIPCHI2,*dsSB);

	RooRealVar fSigInPeak("fSigInPeak","",fsig_1/fsig_0);
	RooRealVar fBkgInPeak("fBkgInPeak","",fbkg_1/fbkg_0);
	RooRealVar fPrompt("fPrompt","frac. prompt",0.5,0.,1.);
	RooFormulaVar promptYield("promptYield","","@0*@1*@2",RooArgList(fPrompt,fSigInPeak,sigYield));
	RooFormulaVar displacedYield("displacedYield","","(1.0-@0)*@1*@2",RooArgList(fPrompt,fSigInPeak,sigYield));
	RooFormulaVar bkgInPeakYield("bkgInPeakYield","","@0*@1",RooArgList(fBkgInPeak,bkgYield));

	RooAddPdf data_pdf2( "data_pdf2",  "data_pdf2", RooArgList(promptLOGIPCHI2,displacedLOGIPCHI2,bkgLOGIPCHI2), RooArgList(promptYield,displacedYield,bkgInPeakYield) );
	RooProdPdf data_pdf2_full("data_pdf2_full", "data_pdf2_full", RooArgList(data_pdf2,bkgYield_constraint));

	if(fix) {
		//promptMean.setConstant();
		//promptWidth.setConstant();
		promptAsym.setConstant();
		promptRhoL.setConstant();
		promptRhoR.setConstant();
		//displacedMean.setConstant();
		//displacedWidth.setConstant();
	}
	if(_sample==charmSample) {
		fPrompt.setVal(1.0);
		fPrompt.setConstant();
		displacedMean.setConstant();
		displacedWidth.setConstant();
	}
	if(_sample==beautySample) {
		fPrompt.setVal(0.0);
		fPrompt.setConstant();
		promptMean.setConstant();
		promptWidth.setConstant();
	}

	gSystem->RedirectOutput(gSaveDir+"/"+dType+"_"+_name+"_"+ptStr+typeName(_type)+"_fits.log","a");
	RooFitResult * result2 = data_pdf2_full.fitTo( *dsPeak, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE), RooFit::Constrain(RooArgSet(bkgYield)));
	gSystem->RedirectOutput(0);

	double fsigIPCHI2_1 = promptLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("sigIPCHI2"))->getVal();
	double fsigIPCHI2_0 = promptLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("fullIPCHI2"))->getVal();

	double f2ndIPCHI2_1 = displacedLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("sigIPCHI2"))->getVal();
	double f2ndIPCHI2_0 = displacedLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("fullIPCHI2"))->getVal();

	double fbkgIPCHI2_1 = bkgLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("sigIPCHI2"))->getVal();
	double fbkgIPCHI2_0 = bkgLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("fullIPCHI2"))->getVal();

	std::cout << "yields in IP window" << std::endl;
	std::cout << promptYield.getVal()*fsigIPCHI2_1/fsigIPCHI2_0 << "\t" << displacedYield.getVal()*f2ndIPCHI2_1/f2ndIPCHI2_0 << "\t" << bkgInPeakYield.getVal()*fbkgIPCHI2_1/fbkgIPCHI2_0 << std::endl;

	//extract results - corrected for the mass window cut
	if(_flavour==5) {
		yield = (fsig_0/fsig_1)*displacedYield.getVal();
		error = (fsig_0/fsig_1)*displacedYield.getPropagatedError(*result2);
	} else {
		yield = (fsig_0/fsig_1)*promptYield.getVal();
		error = (fsig_0/fsig_1)*promptYield.getPropagatedError(*result2);
	}
	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "sigMass" );
	std::vector<std::string> bkg_pdfs;
	bkg_pdfs.push_back( "bkgMass" );

	std::vector<std::string> sig_pdfs2;
	sig_pdfs2.push_back( "displacedLOGIPCHI2" );
	sig_pdfs2.push_back( "promptLOGIPCHI2" );
	std::vector<std::string> bkg_pdfs2;
	bkg_pdfs2.push_back( "bkgLOGIPCHI2" );

	plotFit(DM, massVal-massScaleFact*80., massVal+massScaleFact*80., 40, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+_name+"_M_"+ptStr, "m");
	plotFit(DLOGIPCHI2, -5., 15., 20, dsPeak, data_pdf2, sig_pdfs2, bkg_pdfs2, dType+"_"+_name+"_IPChi2_"+ptStr, "log(IP#chi^{2})");

	////print parameters
	RooArgList params;
	params.add(sigYield);
	params.add(bkgYield);
	params.add(fPrompt);
	params.add(dMass);
	params.add(dWidth);
	params.add(p0);
	params.add(promptMean);
	params.add(promptWidth);
	//params.add(promptAsym);
	//params.add(promptRhoL);
	//params.add(promptRhoR);
	params.add(displacedMean);
	params.add(displacedWidth);
	printParams(dType+"_"+_name+"_params_"+ptStr+".dat",params);

	return true;
}

bool DFitter::fitD_SBS(double& yield, double& error) {
	std::cout << "INFO : fitting D0 - sideband subtraction and IP cut to file " << _name << ", pT in (" << _jptmin << "," << _jptmax << ")"  << std::endl;
	TString ptStr;
	ptStr+=_jptmin; ptStr+="-"; ptStr+=_jptmax;

	if(_dptmin!=5000. || _dptmax!=-1) {
		ptStr+="_D0";
		ptStr+=_dptmin;
		ptStr+="-";
		ptStr+=_dptmax;
	}
	if(!_useEffs) {
		ptStr+="_noEffCorr";
	}

	TFile* f0 = TFile::Open(dFileName());

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TString weightName("");

	if(_flavour==4) {
		if(_useEffs) weightName = "weight4";
	} else if(_flavour==5) {
		if(_useEffs) weightName = "weight5";
	} else {
		return false;
	}

	TString dType("D0");
	double massVal(1865.);
	double massScaleFact(1.);

	//D0 settings below
	massVal=1865.;

	TH1D hpeak("hpeak","",40,-5.,15.);
	TH1D hside("hside","",40,-5.,15.);
	TH1D hsub ("hsub" ,"",40,-5.,15.);

	hpeak.Sumw2();
	hside.Sumw2();
	hsub.Sumw2();

	TString jetPTCut = "(";
	jetPTCut+="JetPT>"; jetPTCut+=_jptmin; jetPTCut+=" && ";
	jetPTCut+="JetPT<"; jetPTCut+=_jptmax; jetPTCut+=")";

	TString peakCut("");
	if(_useEffs) peakCut+=weightName+" * ";
	peakCut+="(TMath::Abs(";
	peakCut+=dType;
	peakCut+="M-";
	peakCut+=massVal;
	peakCut+=")<";
	peakCut+=massScaleFact*30.;
	peakCut+=") * ";
	peakCut+=jetPTCut;
	TString sideCut("");
	if(_useEffs) sideCut+=weightName+" * ";
	sideCut+="(TMath::Abs(";
	sideCut+=dType;
	sideCut+="M-";
	sideCut+=massVal;
	sideCut+=")<";
	sideCut+=massScaleFact*80.;
	sideCut+=" && ";
	sideCut+="TMath::Abs(";
	sideCut+=dType;
	sideCut+="M-";
	sideCut+=massVal;
	sideCut+=")>";
	sideCut+=massScaleFact*50.;
	sideCut+=") * ";
	sideCut+=jetPTCut;

	t0->Draw(dType+"LogIPChi2>>hpeak",peakCut);
	t0->Draw(dType+"LogIPChi2>>hside",sideCut);

	hsub.Add(&hpeak,&hside,1.,-1.);

	TCanvas c;
	hpeak.Draw();
	hsub.SetLineColor(kGreen+2);
	hsub.Draw("same");
	hside.SetLineColor(kRed);
	hside.Draw("same");
	c.SaveAs(gSaveDir+"/"+dType+"_"+_name+"_IPChi2_SBS_"+ptStr+".pdf");

	std::cout << hsub.Integral() << "\t" << hsub.Integral(1,14) << "\t" << hsub.Integral(15,40) << std::endl;

	//extract results - not yet corrected for prompt<->displaced migration
	yield=0.;
	error=0.;
	if(_flavour==5) {
		for(int i=15; i<41; ++i) {
			yield += hsub.GetBinContent(i);
			error += hsub.GetBinError(i)*hsub.GetBinError(i);
		}
		error=TMath::Sqrt(error);
	} else {
		for(int i=1; i<15; ++i) {
			yield += hsub.GetBinContent(i);
			error += hsub.GetBinError(i)*hsub.GetBinError(i);
		}
		error=TMath::Sqrt(error);
	}

	return true;
}

bool DFitter::fitD_SBS1D(double& yield, double& error) {
	std::cout << "INFO : fitting D0 - sideband subtraction and IP fit to file " << _name << ", pT in (" << _jptmin << "," << _jptmax << ")"  << std::endl;
	TString ptStr;
	ptStr+=_jptmin; ptStr+="-"; ptStr+=_jptmax;

	if(_dptmin!=5000. || _dptmax!=-1) {
		ptStr+="_D0";
		ptStr+=_dptmin;
		ptStr+="-";
		ptStr+=_dptmax;
	}
	if(!_useEffs) {
		ptStr+="_noEffCorr";
	}

	TFile* f0 = TFile::Open(dFileName());

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TString weightName="";

	if(_flavour==4) {
		if(_useEffs) weightName = "weight4";
	} else if(_flavour==5) {
		if(_useEffs) weightName = "weight5";
	} else {
		return false;
	}

	TString dType("D0");
	double massVal(1865.);
	double massScaleFact(1.);

	//D0 settings below
	massVal=1865.;

	int nbins(40);

	TH1D hpeak("hpeak","",nbins,-5.,15.);
	TH1D hside("hside","",nbins,-5.,15.);
	TH1D hsub ("hsub" ,"",nbins,-5.,15.);

	TH1D mpeak("mpeak","",80,massVal-massScaleFact*80.,massVal+massScaleFact*80.);
	TH1D mside("mside","",80,massVal-massScaleFact*80.,massVal+massScaleFact*80.);
	TH1D mall ("mall" ,"",80,massVal-massScaleFact*80.,massVal+massScaleFact*80.);

	hpeak.Sumw2();
	hside.Sumw2();
	hsub.Sumw2();

	TString jetPTCut = "(";
	jetPTCut+="JetPT>"; jetPTCut+=_jptmin; jetPTCut+=" && ";
	jetPTCut+="JetPT<"; jetPTCut+=_jptmax; jetPTCut+=")";

	TString peakCut("");
	if(_useEffs) peakCut+=weightName+" * ";
	peakCut+="(TMath::Abs(";
	peakCut+=dType;
	peakCut+="M-";
	peakCut+=massVal;
	peakCut+=")<";
	peakCut+=massScaleFact*40.;
	peakCut+=") * ";
	peakCut+=jetPTCut;
	TString sideCut("");
	if(_useEffs) sideCut+=weightName+" * ";
	sideCut+="(TMath::Abs(";
	sideCut+=dType;
	sideCut+="M-";
	sideCut+=massVal;
	sideCut+=")<";
	sideCut+=massScaleFact*80.;
	sideCut+=" && ";
	sideCut+="TMath::Abs(";
	sideCut+=dType;
	sideCut+="M-";
	sideCut+=massVal;
	sideCut+=")>";
	sideCut+=massScaleFact*40.;
	sideCut+=") * ";
	sideCut+=jetPTCut;
	TString allCut("");
	if(_useEffs) allCut+=weightName+" * ";
	allCut+="(TMath::Abs(";
	allCut+=dType;
	allCut+="M-";
	allCut+=massVal;
	allCut+=")<";
	allCut+=massScaleFact*80.;
	allCut+=") * ";
	allCut+=jetPTCut;

	t0->Draw(dType+"LogIPChi2>>hpeak",peakCut);
	t0->Draw(dType+"LogIPChi2>>hside",sideCut);

	t0->Draw(dType+"M>>mpeak",peakCut);
	t0->Draw(dType+"M>>mside",sideCut);
	t0->Draw(dType+"M>>mall", allCut);

	hsub.Add(&hpeak,&hside,1.,-1.);
	//TODO for(int ibin=1; ibin<=hsub.GetNbinsX(); ++ibin) {
	//TODO 	if(hsub.GetBinContent(ibin)<0.) {
	//TODO 		std::cout << "Resetting negative bin " << ibin << " (" << hsub.GetBinContent(ibin) << " +/- " << hsub.GetBinError(ibin) << ") to 0+/-1." << std::endl;
	//TODO 		hsub.SetBinContent(ibin,0.);
	//TODO 		hsub.SetBinError(ibin,1.);
	//TODO 	}
	//TODO }

	TString axisStr = "Events / ("; axisStr+=massScaleFact*2; axisStr+=" MeV/#it{c}^{2})";
	TCanvas c;
	mpeak.GetXaxis()->SetTitle("#it{m}(#it{D}^{0}) [MeV/#it{c}^{2}]");
	mpeak.GetYaxis()->SetTitle(axisStr);
	mpeak.SetLineColor(kBlue);
	mpeak.SetFillColor(kBlue);
	mpeak.SetFillStyle(3004);
	mpeak.Draw();
	mside.SetLineColor(kRed);
	mside.SetFillColor(kRed);
	mside.SetFillStyle(3005);
	mside.Draw("same");
	mall.SetLineColor(kBlack);
	mall.Draw("hist same");
	c.SaveAs(gSaveDir+"/"+dType+"_"+_name+"_M_SBS_"+ptStr+".pdf");

	hpeak.GetXaxis()->SetTitle("log #it{#chi}^{2}_{IP}(#it{D}^{0})");
	hpeak.GetYaxis()->SetTitle("Events / (0.5)");
	hpeak.SetLineColor(kBlue);
	hpeak.Draw("E1");
	hsub.SetLineColor(kGreen+2);
	hsub.Draw("E1 same");
	hside.SetLineColor(kRed);
	hside.Draw("E1 same");
	c.SaveAs(gSaveDir+"/"+dType+"_"+_name+"_IPChi2_SBS_"+ptStr+".pdf");

	//now fit log(chi^2_IP)

	double valPromptMean(0.), valPromptWidth(0.), valPromptAsym(0.), valPromptRhoL(0.), valPromptRhoR(0.),  valDisplacedMean(5.), valDisplacedWidth(1.5);
	//D0 settings
	valPromptMean = 0.87;
	valPromptWidth = 1.09;
	valPromptAsym = -0.29;
	valPromptRhoL = 1.30;
	valPromptRhoR = 1.69;

	RooRealVar DLOGIPCHI2(  dType+"LogIPChi2",  dType+"LogIPChi2",  -5.,  15.,  ""); 
	DLOGIPCHI2.setBins(nbins);

	RooDataHist dh("dh", "dh", RooArgList(DLOGIPCHI2), &hsub);

	RooRealVar promptMean("promptMean","mean prompt",valPromptMean,-1.,3.);
	RooRealVar promptWidth("promptWidth","width prompt",valPromptWidth,0.5,5.);
	RooRealVar promptAsym("promptAsym","asym. prompt",valPromptAsym,-1.,1.);
	RooRealVar promptRhoL("promptRhoL","exp L prompt",valPromptRhoL,0.5,3.);
	RooRealVar promptRhoR("promptRhoR","exp R prompt",valPromptRhoR,0.5,3.);
	RooPromptShape prompt("prompt","",DLOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	bool fix(true);
	if(fix) {
		//	promptMean.setConstant();
		//	promptWidth.setConstant();
		promptAsym.setConstant();
		promptRhoL.setConstant();
		promptRhoR.setConstant();
		//	displacedMean.setConstant();
		//	displacedWidth.setConstant();
	}

	RooRealVar displacedMean("displacedMean","mean displ.",valDisplacedMean,2.,12.);
	RooRealVar displacedWidth("displacedWidth","width displ.",valDisplacedWidth,0.5,4.);
	RooGaussian displaced("displaced","",DLOGIPCHI2,displacedMean,displacedWidth);

	double maxyield = dh.sumEntries();
	// -- yields
	RooRealVar promptYield(     "promptYield",     "", maxyield/2.,    0.0,  1.1*maxyield);
	RooRealVar displacedYield(  "displacedYield",  "", maxyield/2.,    0.0,  1.1*maxyield);

	if(_sample==charmSample) {
		displacedYield.setVal(0.0);
		displacedYield.setConstant();
		displacedMean.setConstant();
		displacedWidth.setConstant();
	}
	if(_sample==beautySample) {
		promptYield.setVal(0.0);
		promptYield.setConstant();
		promptMean.setConstant();
		promptWidth.setConstant();
	}

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(prompt,displaced), RooArgList(promptYield,displacedYield) );

	gSystem->RedirectOutput(gSaveDir+"/"+dType+"_"+_name+"_"+ptStr+typeName(_type)+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));
	gSystem->RedirectOutput(0);

	//extract results - corrected for the mass window cut
	if(_flavour==5) {
		yield = displacedYield.getVal();
		error = displacedYield.getError();
	} else {
		yield = promptYield.getVal();
		error = promptYield.getError();
	}

	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "displaced" );
	sig_pdfs.push_back( "prompt" );
	std::vector<std::string> bkg_pdfs;

	plotFit(DLOGIPCHI2, -5., 15., 20, &dh, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+_name+"_IPChi2_SBS1D_"+ptStr, "log(IP#chi^{2})");

	////print parameters
	RooArgList params;
	params.add(promptYield);
	params.add(displacedYield);
	params.add(promptMean);
	params.add(promptWidth);
	//params.add(promptAsym);
	//params.add(promptRhoL);
	//params.add(promptRhoR);
	params.add(displacedMean);
	params.add(displacedWidth);
	printParams(dType+"_"+_name+"_params_SBS1D_"+ptStr+".dat",params);

	return true;
}

bool DFitter::fitD_2D(double& yield, double& error) {
	std::cout << "INFO : fitting D0 - 2D mass and IP fit to file " << _name << ", pT in (" << _jptmin << "," << _jptmax << ")"  << std::endl;
	bool fix(true);

	double valPromptMean(0.), valPromptWidth(0.), valPromptAsym(0.), valPromptRhoL(0.), valPromptRhoR(0.),  valDisplacedMean(5.), valDisplacedWidth(1.5);
	//D0 settings
	valPromptMean = 0.87;
	valPromptWidth = 1.09;
	valPromptAsym = -0.29;
	valPromptRhoL = 1.30;
	valPromptRhoR = 1.69;

	TString ptStr;
	ptStr+=_jptmin; ptStr+="-"; ptStr+=_jptmax;

	if(_dptmin!=5000. || _dptmax!=-1) {
		ptStr+="_D0";
		ptStr+=_dptmin;
		ptStr+="-";
		ptStr+=_dptmax;
	}
	if(!_useEffs) {
		ptStr+="_noEffCorr";
	}

	TFile* f0 = TFile::Open(dFileName());

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f1 = TFile::Open("../skim/combHist.root");
	TH2D* combhist = dynamic_cast<TH2D*>(f1->Get("comb"));

	if(!t0 || !combhist) return false;

	TString weightName;

	if(_flavour==4) {
		if(_useEffs) weightName = "weight4";
	} else if(_flavour==5) {
		if(_useEffs) weightName = "weight5";
	} else {
		return false;
	}

	TString dType("D0");
	double massVal(1865.);
	double massScaleFact(1.);
	double ratioVal(3.), fracVal(0.9), alphaVal(3.), nVal(1.);

	//D0 settings below
	massVal=1865.;
	ratioVal=2.85;
	fracVal=0.89;
	alphaVal=3.0;
	nVal=1.0;

	RooRealVar DM(  dType+"M",  dType+"M",  massVal-massScaleFact*80.,  massVal+massScaleFact*80.,  "MeV/#it{c}^{2}"); 
	DM.setBins(30);
	RooRealVar DLOGIPCHI2(  dType+"LogIPChi2",  dType+"LogIPChi2",  -5.,  15.,  ""); 
	DLOGIPCHI2.setBins(20);
	RooRealVar DPT(  dType+"PT",  dType+"PT",  0.,  1000000.); 
	RooRealVar JetPT( "JetPT",  "JetPT",  _jptmin,  _jptmax,  "MeV/#it{c}"); 
	RooRealVar DWeight(  weightName,  weightName,  0.,   1000.); 

	RooArgSet obs;
	obs.add(JetPT);
	obs.add(DM);
	obs.add(DPT);
	obs.add(DLOGIPCHI2);
	if(_useEffs) obs.add(DWeight);

	RooDataSet* ds(0);
	if(_useEffs) ds = new RooDataSet("ds","ds", obs, RooFit::WeightVar(weightName), RooFit::Import(*t0));
	else ds = new RooDataSet("ds","ds", obs, RooFit::Import(*t0));

	//mass PDFs
	RooRealVar dMass("dMass","mass mean",massVal,massVal-massScaleFact*20.,massVal+massScaleFact*20.);
	RooRealVar dWidth("dWidth","mass width",massScaleFact*20.,massScaleFact*2.,massScaleFact*100.);

	RooRealVar dRatio("dRatio","ratio of widths", ratioVal, 1.0, 5);
	RooFormulaVar dWidthG("dWidthG","","@0*@1",RooArgList(dWidth,dRatio));
	RooRealVar dAlpha("dAlpha", "dAlpha", alphaVal, 0.5, 5.0);
	RooRealVar dN(    "dN",     "dN",     nVal);//2.0, 0.0, 5.0);
	RooRealVar dFrac("dFrac","fraction in CB", fracVal, 0., 1.);

	RooGaussian sigMass_gauss("sigMass_gauss","",DM,dMass,dWidthG);
	RooCBShape  sigMass_cb("sigMass_cb", "", DM, dMass, dWidth, dAlpha, dN); 
	RooAddPdf sigMass( "sigMass", "", RooArgList(sigMass_cb,sigMass_gauss), RooArgList(dFrac) );

	//IPchi2 PDFs
	RooRealVar promptMean("promptMean","mean prompt",valPromptMean,-1.,3.);
	RooRealVar promptWidth("promptWidth","width prompt",valPromptWidth,0.5,5.);
	RooRealVar promptAsym("promptAsym","asym. prompt",valPromptAsym,-1.,1.);
	RooRealVar promptRhoL("promptRhoL","exp L prompt",valPromptRhoL,0.5,3.);
	RooRealVar promptRhoR("promptRhoR","exp R prompt",valPromptRhoR,0.5,3.);
	RooPromptShape promptLOGIPCHI2("promptLOGIPCHI2","",DLOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	RooRealVar displacedMean("displacedMean","mean displ.",valDisplacedMean,2.,12.);
	RooRealVar displacedWidth("displacedWidth","width displ.",valDisplacedWidth,0.5,4.);
	RooGaussian displacedLOGIPCHI2("displacedLOGIPCHI2","",DLOGIPCHI2,displacedMean,displacedWidth);

	//2D PDFs
	RooDataHist bkgData ("bkgData",   "bkgData",   RooArgList(DM,DLOGIPCHI2), combhist);

	RooProdPdf prompt   ("prompt",    "prompt",    RooArgList(sigMass,promptLOGIPCHI2));
	RooProdPdf displaced("displaced", "displaced", RooArgList(sigMass,displacedLOGIPCHI2));
	RooHistPdf bkg      ("bkg",       "bkg",       RooArgSet(DM,DLOGIPCHI2), bkgData); 

	//constants
	dRatio.setConstant();
	dFrac.setConstant();
	dAlpha.setConstant();
	dN.setConstant();
	if(fix) {
		promptAsym.setConstant();
		promptRhoL.setConstant();
		promptRhoR.setConstant();
	}

	double maxyield = 1e4;
	// -- yields
	RooRealVar promptYield(     "promptYield",     "", 750,    0.0,  maxyield);
	RooRealVar displacedYield(  "displacedYield",  "", 750,    0.0,  maxyield);
	RooRealVar bkgYield(        "bkgYield",        "", 750,    0.0,  maxyield);

	if(_sample==charmSample) {
		displacedYield.setVal(0.0);
		displacedYield.setConstant();
		displacedMean.setConstant();
		displacedWidth.setConstant();
	}
	if(_sample==beautySample) {
		promptYield.setVal(0.0);
		promptYield.setConstant();
		promptMean.setConstant();
		promptWidth.setConstant();
	}

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(prompt,displaced,bkg), RooArgList(promptYield,displacedYield,bkgYield) );

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(gSaveDir+"/"+dType+"_"+_name+"_"+ptStr+typeName(_type)+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));
	gSystem->RedirectOutput(0);

	DM.setRange("signal",massVal-massScaleFact*20.,massVal+massScaleFact*20.);
	DM.setRange("sideLo",massVal-massScaleFact*80.,massVal-massScaleFact*60.);
	DM.setRange("sideHi",massVal+massScaleFact*60.,massVal+massScaleFact*80.);
	DM.setRange("full",  massVal-massScaleFact*80.,massVal+massScaleFact*80.);

	DLOGIPCHI2.setRange("sigIPCHI2",-5.,3.);
	DLOGIPCHI2.setRange("fullIPCHI2",-5.,15.);

	//extract results - corrected for the mass window cut
	if(_flavour==5) {
		yield = displacedYield.getVal();
		error = displacedYield.getError();
	} else {
		yield = promptYield.getVal();
		error = promptYield.getError();
	}
	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "prompt" );
	sig_pdfs.push_back( "displaced" );
	std::vector<std::string> bkg_pdfs;
	bkg_pdfs.push_back( "bkg" );

	plotFit(DM, massVal-massScaleFact*80., massVal+massScaleFact*80., 40, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+_name+"_M_2D_"+ptStr, "m");
	plotFit(DLOGIPCHI2, -5., 15., 20, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+_name+"_IPChi2_2D_"+ptStr, "log(IP#chi^{2})");

	////print parameters
	RooArgList params;
	params.add(promptYield);
	params.add(displacedYield);
	params.add(bkgYield);
	params.add(dMass);
	params.add(dWidth);
	params.add(promptMean);
	params.add(promptWidth);
	//params.add(promptAsym);
	//params.add(promptRhoL);
	//params.add(promptRhoR);
	params.add(displacedMean);
	params.add(displacedWidth);
	printParams(dType+"_"+_name+"_params2D_"+ptStr+".dat",params);

	return true;
}

bool DFitter::addEffs() {
	if(!_inputsSet) {
		std::cout << "ERROR in DFitter::addEffs: input files not set yet" << std::endl;
		return false;
	}

	if(_useSimpleEffs) {
		return addEffsSimple();
	} else {
		return addEffsFull();
	}
}

bool DFitter::testEffs(int flavour) {
	if(!_effsSet) {
		std::cout << "INFO in DFitter::testEffs: D efficiencies not calculated yet" << std::endl;
		std::cout << "                       running now" << std::endl;
		if(!addEffs()) return false;
	}

	_flavour = flavour;

	if(_useSimpleEffs) {
		return testEffsSimple();
	} else {
		return testEffsFull();
	}
}

bool DFitter::fitD(int flavour, double& yield, double& error, uint binPT, uint binY, bool useEffs) {
	if(useEffs && !_effsSet) {
		std::cout << "INFO in DFitter::fitD: D efficiencies not calculated yet" << std::endl;
		std::cout << "                       running now" << std::endl;
		if(!addEffs()) return false;
	}

	_flavour = flavour;
	_jptmin = _ptBins->GetBinLowEdge(binPT);
	_jptmax = _ptBins->GetBinLowEdge(binPT+1);
	if(_yBins) {
		_ymin = _yBins->GetBinLowEdge(binY);
		_ymax = _yBins->GetBinLowEdge(binY+1);
	}
	_useEffs = useEffs;

	switch(_type) {
		case fit1D1D:
			return fitD_1D1D(yield, error);
		case fit2D:
			return fitD_2D(yield, error);
		case fitSidebandSub:
			return fitD_SBS(yield, error);
		case fitSBS1D:
			return fitD_SBS1D(yield, error);
		default:
			return false;
	}
}

bool DFitter::testEffsFull() {
	if(!_dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << _name << std::endl;
	TFile* f0 = TFile::Open(_dataFile);

	TFile* f2 = TFile::Open(_effFile);
	TFile* f3 = TFile::Open(_accFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
	TH2D* hcor = dynamic_cast<TH2D*>(f3->Get("corr"));

	TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	//TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pidmceffD045"));//TODO
	TH2D* hsel(0);
	if(_flavour==5) hsel = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
	else hsel = dynamic_cast<TH2D*>(f2->Get("seleffD04"));
	//TH2D* hacc  = dynamic_cast<TH2D*>(f2->Get("acceffD045"));
	//TH2D* hrec  = dynamic_cast<TH2D*>(f2->Get("receffD045"));

	if(!hacc || !hrec || !hcor || !hpid || !hsel) return false;

	//yield histograms - 0=all; 1=passAcc; 2=passRec; 3=passSel; 4=passPID.
	//tru=truth D0; fit=eff-corrected reco D0; fnd=eff-corrected truth-matched reco D0.
	std::vector<TString> stages;
	stages.push_back("all");
	stages.push_back("geometric");
	stages.push_back("reconstruction");
	stages.push_back("selection");
	stages.push_back("PID");
	std::vector<TH1D*> truD0truePT;//the true distribution at this stage
	truD0truePT.push_back(cloneTH1D("tru0truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru1truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru2truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru3truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru4truept", _ptBins));
	std::vector<TH1D*> effD0truePT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effD0truePT.push_back(cloneTH1D("eff0truept", _ptBins));
	effD0truePT.push_back(cloneTH1D("eff1truept", _ptBins));
	effD0truePT.push_back(cloneTH1D("eff2truept", _ptBins));
	effD0truePT.push_back(cloneTH1D("eff3truept", _ptBins));
	effD0truePT.push_back(cloneTH1D("eff4truept", _ptBins));
	std::vector<TH1D*> fitD0truePT;//sideband-subtracted distribution at this stage
	fitD0truePT.push_back(cloneTH1D("fit0truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit1truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit2truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit3truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit4truept", _ptBins));
	std::vector<TH1D*> fndD0truePT;//truth-matched distribution at this stage
	fndD0truePT.push_back(cloneTH1D("fnd0truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd1truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd2truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd3truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd4truept", _ptBins));
	std::vector<TH1D*> truD0recoPT;
	truD0recoPT.push_back(cloneTH1D("tru0recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru1recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru2recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru3recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru4recopt", _ptBins));
	std::vector<TH1D*> effD0recoPT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effD0recoPT.push_back(cloneTH1D("eff0recopt", _ptBins));
	effD0recoPT.push_back(cloneTH1D("eff1recopt", _ptBins));
	effD0recoPT.push_back(cloneTH1D("eff2recopt", _ptBins));
	effD0recoPT.push_back(cloneTH1D("eff3recopt", _ptBins));
	effD0recoPT.push_back(cloneTH1D("eff4recopt", _ptBins));
	std::vector<TH1D*> fitD0recoPT;
	fitD0recoPT.push_back(cloneTH1D("fit0recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit1recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit2recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit3recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit4recopt", _ptBins));
	std::vector<TH1D*> fndD0recoPT;
	fndD0recoPT.push_back(cloneTH1D("fnd0recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd1recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd2recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd3recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd4recopt", _ptBins));

	//entries in this vector are true and reco for each jet pT bin in ascending order
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

	for(int i=0; i<5; ++i) {
		truD0truePT.at(i)->Sumw2();
		effD0truePT.at(i)->Sumw2();
		fitD0truePT.at(i)->Sumw2();
		fndD0truePT.at(i)->Sumw2();
		truD0recoPT.at(i)->Sumw2();
		effD0recoPT.at(i)->Sumw2();
		fitD0recoPT.at(i)->Sumw2();
		fndD0recoPT.at(i)->Sumw2();
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

	t0->SetBranchAddress("D0M",           &D0M);
	t0->SetBranchAddress("D0IPCHI2",      &D0IPCHI2);
	t0->SetBranchAddress("D0PT",          &D0PT);
	t0->SetBranchAddress("D0PX",          &D0PX);
	t0->SetBranchAddress("D0PY",          &D0PY);
	t0->SetBranchAddress("D0PZ",          &D0PZ);
	t0->SetBranchAddress("D0E",           &D0E);
	t0->SetBranchAddress("D0X",           &D0X);
	t0->SetBranchAddress("D0Y",           &D0Y);
	t0->SetBranchAddress("D0Z",           &D0Z);
	t0->SetBranchAddress("D0KP",          &D0KP);
	t0->SetBranchAddress("D0KPT",         &D0KPT);
	t0->SetBranchAddress("D0KPX",         &D0KPX);
	t0->SetBranchAddress("D0KPY",         &D0KPY);
	t0->SetBranchAddress("D0KPZ",         &D0KPZ);
	t0->SetBranchAddress("D0PIP",         &D0PIP);
	t0->SetBranchAddress("D0PIPT",        &D0PIPT);
	t0->SetBranchAddress("D0PIPX",        &D0PIPX);
	t0->SetBranchAddress("D0PIPY",        &D0PIPY);
	t0->SetBranchAddress("D0PIPZ",        &D0PIPZ);
	t0->SetBranchAddress("D0KWEIGHT",     &D0KWEIGHT);
	t0->SetBranchAddress("D0PIWEIGHT",    &D0PIWEIGHT);
	t0->SetBranchAddress("D0TRUEIDX",     &D0TRUEIDX);
	t0->SetBranchAddress("D0TRUETRK0",    &D0TRUETRK0);
	t0->SetBranchAddress("D0TRUETRK1",    &D0TRUETRK1);

	t0->SetBranchAddress("TRUEDID",       &TRUEDID);
	t0->SetBranchAddress("TRUEDPX",       &TRUEDPX);
	t0->SetBranchAddress("TRUEDPY",       &TRUEDPY);
	t0->SetBranchAddress("TRUEDPZ",       &TRUEDPZ);
	t0->SetBranchAddress("TRUEDE",        &TRUEDE);
	t0->SetBranchAddress("TRUEDX",        &TRUEDX);
	t0->SetBranchAddress("TRUEDY",        &TRUEDY);
	t0->SetBranchAddress("TRUEDZ",        &TRUEDZ);
	t0->SetBranchAddress("TRUEDFROMB",    &TRUEDFROMB);
	t0->SetBranchAddress("TRUEDTRK0P",    &TRUEDTRK0P);
	t0->SetBranchAddress("TRUEDTRK0PT",   &TRUEDTRK0PT);
	t0->SetBranchAddress("TRUEDTRK0INACC",&TRUEDTRK0INACC);
	t0->SetBranchAddress("TRUEDTRK0RECO", &TRUEDTRK0RECO);
	t0->SetBranchAddress("TRUEDTRK1P",    &TRUEDTRK1P);
	t0->SetBranchAddress("TRUEDTRK1PT",   &TRUEDTRK1PT);
	t0->SetBranchAddress("TRUEDTRK1INACC",&TRUEDTRK1INACC);
	t0->SetBranchAddress("TRUEDTRK1RECO", &TRUEDTRK1RECO);
	t0->SetBranchAddress("TRUEDTRK0IDX",  &TRUEDTRK0IDX);
	t0->SetBranchAddress("TRUEDTRK1IDX",  &TRUEDTRK1IDX);
	t0->SetBranchAddress("TRUEDTRK2IDX",  &TRUEDTRK2IDX);
	t0->SetBranchAddress("TRUEDTRK3IDX",  &TRUEDTRK3IDX);
	t0->SetBranchAddress("TRUEDTRK0ID",  &TRUEDTRK0ID);
	t0->SetBranchAddress("TRUEDTRK1ID",  &TRUEDTRK1ID);
	t0->SetBranchAddress("TRUEDTRK2ID",  &TRUEDTRK2ID);
	t0->SetBranchAddress("TRUEDTRK3ID",  &TRUEDTRK3ID);

	t0->SetBranchAddress("JetPT",         &JetPT);
	t0->SetBranchAddress("JetEta",        &JetEta);
	t0->SetBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = t0->GetEntries();

	double dpt(0.), deta(0.);

	boost::progress_display progress( nentries0 );
	for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
		++progress;
		t0->GetEntry(ientry);

		std::vector<int> trueD0s, foundD0s;

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
			if(_useRhoZCorr) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));

			truD0recoPT.at(0)->Fill(JetPT);
			truD0truePT.at(0)->Fill(JetTruePT);
			d0Kinematics[0]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;

			truD0recoPT.at(1)->Fill(JetPT);
			truD0truePT.at(1)->Fill(JetTruePT);
			if(effacc>0.) {
				effD0recoPT.at(0)->Fill(JetPT,1./effacc);
				effD0truePT.at(0)->Fill(JetTruePT,1./effacc);
			} else {
				std::cout << "ACC=" << effacc << " at " << dpt << "," << deta << std::endl;
			}
			d0Kinematics[2]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			truD0recoPT.at(2)->Fill(JetPT);
			truD0truePT.at(2)->Fill(JetTruePT);
			if(effrec>0.) {
				effD0recoPT.at(1)->Fill(JetPT,1./effrec);
				effD0truePT.at(1)->Fill(JetTruePT,1./effrec);
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
					if(effsel>0.) {
						effD0recoPT.at(2)->Fill(JetPT,1./effsel);
						effD0truePT.at(2)->Fill(JetTruePT,1./effsel);
					} else {
						std::cout << "SEL=" << effsel << " at " << dpt << "," << deta << std::endl;
					}
					d0Kinematics[6]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

					//use weights for PID
					double weight = D0KWEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
					if(D0P0.Pt()>25000. || D0P0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT
					truD0recoPT.at(4)->Fill(JetPT,weight);
					truD0truePT.at(4)->Fill(JetTruePT,weight);
					if(effpid>0.) {
						effD0recoPT.at(3)->Fill(JetPT,weight/effpid);
						effD0truePT.at(3)->Fill(JetTruePT,weight/effpid);
					} else {
						std::cout << "PID=" << effpid << " at " << dpt << "," << deta << std::endl;
					}
					d0Kinematics[8]->Fill(TRUEDP.Pt(),TRUEDP.Eta(),weight);

					trueD0s.push_back(d*100+s);

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
			if(_useRhoZCorr) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));
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
			if(!matched) std::cout << ientry << ": truD0 " << trueD0s[i] << " has no matching fndD0" << std::endl;
		}
		for(unsigned int i=0; i<foundD0s.size(); ++i) {
			bool matched(false);
			for(unsigned int j=0; j<trueD0s.size(); ++j) {
				if(foundD0s[i]==trueD0s[j]) {
					matched=true;
					break;
				}
			}
			if(!matched) std::cout << ientry << ": fndD0 " << trueD0s[i] << " has no matching truD0" << std::endl;
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

	return true;
}

bool DFitter::testEffsSimple() {
	if(!_dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << _name << std::endl;
	TFile* f0 = TFile::Open(_dataFile);

	TFile* f2 = TFile::Open(_effFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* heff(0);
	if(_flavour==5) heff = dynamic_cast<TH2D*>(f2->Get("effD05"));
	else heff = dynamic_cast<TH2D*>(f2->Get("effD04"));

	if(!heff) return false;

	//yield histograms - 0=all; 1=passAcc; 2=passRec; 3=passSel; 4=passPID.
	//tru=truth D0; fit=eff-corrected reco D0; fnd=eff-corrected truth-matched reco D0.
	std::vector<TString> stages;
	stages.push_back("all");
	stages.push_back("geometric");
	stages.push_back("reconstruction");
	stages.push_back("selection");
	stages.push_back("PID");
	std::vector<TH1D*> truD0truePT;
	truD0truePT.push_back(cloneTH1D("tru0truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru1truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru2truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru3truept", _ptBins));
	truD0truePT.push_back(cloneTH1D("tru4truept", _ptBins));
	std::vector<TH1D*> fitD0truePT;
	fitD0truePT.push_back(cloneTH1D("fit0truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit1truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit2truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit3truept", _ptBins));
	fitD0truePT.push_back(cloneTH1D("fit4truept", _ptBins));
	std::vector<TH1D*> fndD0truePT;
	fndD0truePT.push_back(cloneTH1D("fnd0truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd1truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd2truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd3truept", _ptBins));
	fndD0truePT.push_back(cloneTH1D("fnd4truept", _ptBins));
	std::vector<TH1D*> truD0recoPT;
	truD0recoPT.push_back(cloneTH1D("tru0recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru1recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru2recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru3recopt", _ptBins));
	truD0recoPT.push_back(cloneTH1D("tru4recopt", _ptBins));
	std::vector<TH1D*> fitD0recoPT;
	fitD0recoPT.push_back(cloneTH1D("fit0recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit1recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit2recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit3recopt", _ptBins));
	fitD0recoPT.push_back(cloneTH1D("fit4recopt", _ptBins));
	std::vector<TH1D*> fndD0recoPT;
	fndD0recoPT.push_back(cloneTH1D("fnd0recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd1recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd2recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd3recopt", _ptBins));
	fndD0recoPT.push_back(cloneTH1D("fnd4recopt", _ptBins));

	//entries in this vector are true and reco for each jet pT bin in ascending order
	std::vector<TH2D*> d0Kinematics;
	d0Kinematics.push_back(cloneTH2D("d0KineTrue0", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineReco0", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue1", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineReco1", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue2", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineReco2", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue3", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineReco3", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineTrue4", heff));
	d0Kinematics.push_back(cloneTH2D("d0KineReco4", heff));
	for(unsigned int i=0; i<d0Kinematics.size(); ++i) {
		d0Kinematics[i]->Sumw2();
	}
	std::vector<TH1D*> d0KinematicsPulls;
	d0KinematicsPulls.push_back(new TH1D("pulls1", "", 30, -1., 1.));
	d0KinematicsPulls.push_back(new TH1D("pulls2", "", 30, -1., 1.));
	d0KinematicsPulls.push_back(new TH1D("pulls3", "", 30, -1., 1.));
	d0KinematicsPulls.push_back(new TH1D("pulls4", "", 30, -1., 1.));

	for(int i=0; i<5; ++i) {
		truD0truePT.at(i)->Sumw2();
		fitD0truePT.at(i)->Sumw2();
		fndD0truePT.at(i)->Sumw2();
		truD0recoPT.at(i)->Sumw2();
		fitD0recoPT.at(i)->Sumw2();
		fndD0recoPT.at(i)->Sumw2();
	}

	std::vector<TH1D*> d0masstru;
	std::vector<TH1D*> d0massfit;
	std::vector<TH1D*> d0massfnd;
	for(int i=0; i< _ptBins->GetNbinsX(); ++i) {
		TString suffix("_");
		suffix+=i;
		d0masstru.push_back(new TH1D("d0masstru"+suffix, "", 80, 1865.-80., 1865.+80.));
		d0massfit.push_back(new TH1D("d0massfit"+suffix, "", 80, 1865.-80., 1865.+80.));
		d0massfnd.push_back(new TH1D("d0massfnd"+suffix, "", 80, 1865.-80., 1865.+80.));
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

	t0->SetBranchAddress("D0M",           &D0M);
	t0->SetBranchAddress("D0IPCHI2",      &D0IPCHI2);
	t0->SetBranchAddress("D0PT",          &D0PT);
	t0->SetBranchAddress("D0PX",          &D0PX);
	t0->SetBranchAddress("D0PY",          &D0PY);
	t0->SetBranchAddress("D0PZ",          &D0PZ);
	t0->SetBranchAddress("D0E",           &D0E);
	t0->SetBranchAddress("D0X",           &D0X);
	t0->SetBranchAddress("D0Y",           &D0Y);
	t0->SetBranchAddress("D0Z",           &D0Z);
	t0->SetBranchAddress("D0KP",          &D0KP);
	t0->SetBranchAddress("D0KPT",         &D0KPT);
	t0->SetBranchAddress("D0KPX",         &D0KPX);
	t0->SetBranchAddress("D0KPY",         &D0KPY);
	t0->SetBranchAddress("D0KPZ",         &D0KPZ);
	t0->SetBranchAddress("D0PIP",         &D0PIP);
	t0->SetBranchAddress("D0PIPT",        &D0PIPT);
	t0->SetBranchAddress("D0PIPX",        &D0PIPX);
	t0->SetBranchAddress("D0PIPY",        &D0PIPY);
	t0->SetBranchAddress("D0PIPZ",        &D0PIPZ);
	t0->SetBranchAddress("D0KWEIGHT",     &D0KWEIGHT);
	t0->SetBranchAddress("D0PIWEIGHT",    &D0PIWEIGHT);
	t0->SetBranchAddress("D0TRUEIDX",     &D0TRUEIDX);
	t0->SetBranchAddress("D0TRUETRK0",    &D0TRUETRK0);
	t0->SetBranchAddress("D0TRUETRK1",    &D0TRUETRK1);

	t0->SetBranchAddress("TRUEDID",       &TRUEDID);
	t0->SetBranchAddress("TRUEDPX",       &TRUEDPX);
	t0->SetBranchAddress("TRUEDPY",       &TRUEDPY);
	t0->SetBranchAddress("TRUEDPZ",       &TRUEDPZ);
	t0->SetBranchAddress("TRUEDE",        &TRUEDE);
	t0->SetBranchAddress("TRUEDFROMB",    &TRUEDFROMB);
	t0->SetBranchAddress("TRUEDTRK0P",    &TRUEDTRK0P);
	t0->SetBranchAddress("TRUEDTRK0PT",   &TRUEDTRK0PT);
	t0->SetBranchAddress("TRUEDTRK0INACC",&TRUEDTRK0INACC);
	t0->SetBranchAddress("TRUEDTRK0RECO", &TRUEDTRK0RECO);
	t0->SetBranchAddress("TRUEDTRK1P",    &TRUEDTRK1P);
	t0->SetBranchAddress("TRUEDTRK1PT",   &TRUEDTRK1PT);
	t0->SetBranchAddress("TRUEDTRK1INACC",&TRUEDTRK1INACC);
	t0->SetBranchAddress("TRUEDTRK1RECO", &TRUEDTRK1RECO);
	t0->SetBranchAddress("TRUEDTRK0IDX",  &TRUEDTRK0IDX);
	t0->SetBranchAddress("TRUEDTRK1IDX",  &TRUEDTRK1IDX);
	t0->SetBranchAddress("TRUEDTRK2IDX",  &TRUEDTRK2IDX);
	t0->SetBranchAddress("TRUEDTRK3IDX",  &TRUEDTRK3IDX);
	t0->SetBranchAddress("TRUEDTRK0ID",  &TRUEDTRK0ID);
	t0->SetBranchAddress("TRUEDTRK1ID",  &TRUEDTRK1ID);
	t0->SetBranchAddress("TRUEDTRK2ID",  &TRUEDTRK2ID);
	t0->SetBranchAddress("TRUEDTRK3ID",  &TRUEDTRK3ID);

	t0->SetBranchAddress("JetPT",         &JetPT);
	t0->SetBranchAddress("JetEta",        &JetEta);
	t0->SetBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = t0->GetEntries();

	double dpt(0.), deta(0.);

	boost::progress_display progress( nentries0 );
	for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
		++progress;
		t0->GetEntry(ientry);

		int ptBin = _ptBins->FindBin(JetPT);
		if(ptBin<1||ptBin>_ptBins->GetNbinsX()) continue;

		for(unsigned int d=0; d<TRUEDID->size(); ++d) {
			//only use D0->Kpi with Pt>5GeV
			if(TMath::Abs(TRUEDID->at(d))!=421) continue;
			if(TRUEDTRK0IDX->at(d)==-1) continue;
			if(TRUEDTRK2IDX->at(d)!=-1) continue;
			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;

			TVector3 TRUEDP (TRUEDPX->at(d)  ,TRUEDPY->at(d)  ,TRUEDPZ->at(d));

			truD0recoPT.at(0)->Fill(JetPT);
			truD0truePT.at(0)->Fill(JetTruePT);
			d0Kinematics[0]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;

			truD0recoPT.at(1)->Fill(JetPT);
			truD0truePT.at(1)->Fill(JetTruePT);
			d0Kinematics[2]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			truD0recoPT.at(2)->Fill(JetPT);
			truD0truePT.at(2)->Fill(JetTruePT);
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
					d0Kinematics[6]->Fill(TRUEDP.Pt(),TRUEDP.Eta());

					//use weights for PID
					double weight = D0KWEIGHT->at(s);// pion PID removed *D0PIWEIGHT->at(s);
					if(D0P0.Pt()>25000. || D0P0.Mag()>500000.) weight = 1.; //PID turned off for high P or PT
					truD0recoPT.at(4)->Fill(JetPT,weight);
					truD0truePT.at(4)->Fill(JetTruePT,weight);
					d0Kinematics[8]->Fill(TRUEDP.Pt(),TRUEDP.Eta(),weight);
					d0masstru.at(ptBin-1)->Fill(D0M->at(s),weight);

					break;
				}
			}
		}
		for(unsigned int s=0; s<D0M->size(); ++s) {
			//first check in our tight acceptance
			TVector3 D0P (D0PX->at(s)  ,D0PY->at(s)  ,D0PZ->at(s));
			TVector3 D0P0(D0KPX->at(s) ,D0KPY->at(s) ,D0KPZ->at(s));
			TVector3 D0P1(D0PIPX->at(s),D0PIPY->at(s),D0PIPZ->at(s));

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
			double eff(0.);

			if(dpt>=100000.) dpt=99999.;
			eff    =      heff ->GetBinContent(heff ->FindBin(dpt      ,deta));

			if(eff<0.01) {
				std::cout << dpt << "\t" << deta << "\t" << eff << std::endl;
				continue;
			}

			double sbSubSwitch(0.);
			if(TMath::Abs(D0M->at(s)-1865)<40.) sbSubSwitch = 1.;
			else if(TMath::Abs(D0M->at(s)-1865)<80.) sbSubSwitch = -1.;

			fitD0recoPT.at(4)->Fill(JetPT,weight*sbSubSwitch);
			fitD0recoPT.at(0)->Fill(JetPT,weight*sbSubSwitch/(eff));

			fitD0truePT.at(4)->Fill(JetTruePT,weight*sbSubSwitch);
			fitD0truePT.at(0)->Fill(JetTruePT,weight*sbSubSwitch/(eff));

			d0Kinematics[9]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch);
			d0Kinematics[1]->Fill(D0P.Pt(),D0P.Eta(),weight*sbSubSwitch/(eff));
			d0massfit.at(ptBin-1)->Fill(D0M->at(s),weight);

			//truth match to real, accepted and reconstructed D0
			if(D0TRUEIDX->at(s)<0) break;
			int d = D0TRUEIDX->at(s);
			if(TRUEDTRK0IDX->at(d)==-1) continue;
			if(TRUEDTRK2IDX->at(d)!=-1) continue;
			if(TRUEDPX->at(d)*TRUEDPX->at(d)+TRUEDPY->at(d)*TRUEDPY->at(d)<5000.*5000.) continue;
			if(TRUEDTRK0INACC->at(d)!=1 || TRUEDTRK1INACC->at(d)!=1) continue;
			if(TRUEDTRK0RECO->at(d)!=1 || TRUEDTRK1RECO->at(d)!=1) continue;

			fndD0recoPT.at(4)->Fill(JetPT,weight);
			fndD0recoPT.at(0)->Fill(JetPT,weight/(eff));

			fndD0truePT.at(4)->Fill(JetTruePT,weight);
			fndD0truePT.at(0)->Fill(JetTruePT,weight/(eff));

			d0massfnd.at(ptBin-1)->Fill(D0M->at(s),weight);

			break;//only keep one D0 candidate per entry
		}
	}

	printf("bins of true jet pT\n");
	for(int i=0; i<5; ++i) {
		printf("stage %d (%s)\n",i,stages[i].Data());
		for(int j=1; j<=truD0truePT.at(i)->GetNbinsX(); ++j) {
			printf("%6.0f-%6.0f\t",truD0truePT.at(i)->GetBinLowEdge(j),truD0truePT.at(i)->GetBinLowEdge(j+1));
			printf("%5.0f\t%5.0f\t%5.0f",truD0truePT.at(i)->GetBinContent(j),fndD0truePT.at(i)->GetBinContent(j),fitD0truePT.at(i)->GetBinContent(j));
			if(i>0) {
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
			printf("%5.0f\t%5.0f\t%5.0f\t",truD0recoPT.at(i)->GetBinContent(j),fndD0recoPT.at(i)->GetBinContent(j),fitD0recoPT.at(i)->GetBinContent(j));
			if(i>0) {
				printf("\t%5.3f\t%5.3f", truD0recoPT.at(i)->GetBinContent(j)/truD0recoPT.at(i-1)->GetBinContent(j),
						         fndD0recoPT.at(i)->GetBinContent(j)/fndD0recoPT.at(i-1)->GetBinContent(j));
				printf("\t%5.3f", (fndD0recoPT.at(i)->GetBinContent(j)/fndD0recoPT.at(i-1)->GetBinContent(j))/
						  (truD0recoPT.at(i)->GetBinContent(j)/truD0recoPT.at(i-1)->GetBinContent(j)));
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
	d0Kinematics[9]->Divide(d0Kinematics[1]);
	d0Kinematics[8]->Divide(d0Kinematics[0]);
	d0Kinematics[9]->Divide(d0Kinematics[8]);
	d0Kinematics[9]->SetMinimum(0.);
	d0Kinematics[9]->SetMaximum(2.);
	d0Kinematics[9]->Draw("colz");
	c1.SaveAs(gSaveDir+"/D0KineRatio.pdf");
	c1.SetLogx(0);
	
	for(int i=0; i< _ptBins->GetNbinsX(); ++i) {
		TString suffix("_");
		suffix+=i;
		double sbSub = (d0massfit.at(i)->Integral(1,20) + d0massfit.at(i)->Integral(61,80))/40.;
		for(int j=1; j<81; ++j) {
			d0massfit.at(i)->SetBinContent(j, d0massfit.at(i)->GetBinContent(j) - sbSub);
		}
		std::cout << d0masstru.at(i)->Integral() << "\t" << d0massfit.at(i)->Integral() << "\t" << d0massfnd.at(i)->Integral() << std::endl;
		std::cout << d0masstru.at(i)->Integral(21,60) << "\t" << d0massfit.at(i)->Integral(21,60) << "\t" << d0massfnd.at(i)->Integral(21,60) << std::endl;
		d0masstru[i]->SetMaximum(1.1*TMath::Max(TMath::Max(d0masstru[i]->GetMaximum(),d0massfit[i]->GetMaximum()),d0massfnd[i]->GetMaximum()));
		d0massfit[i]->SetLineColor(kRed);
		d0massfnd[i]->SetLineColor(kMagenta);
		d0masstru[i]->Draw();
		d0massfit[i]->Draw("same");
		d0massfnd[i]->Draw("same");
		c1.SaveAs(gSaveDir+"/D0Mass"+suffix+".pdf");
	}
	
	return true;
}

TString const DFitter::dFileName() {
	return gSaveDir+"/D0jets"+_name+".root";
}

TString const DFitter::typeName(fitType type) {
	switch(type) {
		case fit1D1D:
			return "1D1Dfit";
		case fit2D:
			return "2Dfit";
		case fitSidebandSub:
			return "SBSSBSfit";
		case fitSBS1D:
			return "SBS1Dfit";
		default:
			return "";
	}
}
