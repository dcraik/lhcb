#include <vector>
#include <fstream>

#include "TCanvas.h"
#include "TChain.h"
#include "TError.h"
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
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include "RooStats/SPlot.h"

#include "RooPromptShape.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"

#include <boost/progress.hpp>

enum fitType{
	fit1D1D,
	fit2D,
	fitSidebandSub,
	fitSBS1D
};

enum svFitType{
	fitBoth,
	fitMCorr,
	fitNTrk
};

//globals to save passing these around
fitType whichFit(fitSBS1D);
svFitType whichSVFit(fitBoth);
double d0minpt(5000);
double d0maxpt(-1);
TString savedir("output");
bool useSimpleEff(false);
bool useRhoZEffCor(true);

//the following globals give the locations of input tuples
//these may be overridden in certain cases, e.g. if doing an MC closure test
TString charmSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_14X.root";
TString beautySimFile = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_15X.root";
TString lightSimFile  = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_101.root";
TString dataFile      = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_100.root";

//efficiency inputs - update these with latest version numbers
TString simpleEffFile = "../efficiencies/SimpleEffs_2XX_16x8bins_up190216.root";
TString effFile = "../efficiencies/D0Effs_2XX_16x8bins_up190213.root";
TString accFile = "../efficiencies/D0AccEffNewUp190205.root";

bool dataIsMC(false);
bool dataIsResampledMC(false);

enum jetType {
	jetRecoD04,
	jetRecoD05,
	jetRecoSV4,
	jetRecoSV5
};

//functions to clone a histogram
TH1D* cloneTH1D(TString name, TH1D* hist)
{
	if (!hist) return 0;
	TH1D* cloneHist = new TH1D(name,
			"",
			hist->GetNbinsX(),
			hist->GetXaxis()->GetXbins()->GetArray());
	return cloneHist;
}

TH2D* cloneTH2D(TString name, TH2D* hist)
{
	if (!hist) return 0;
	TH2D* cloneHist = new TH2D(name,
			"",
			hist->GetNbinsX(),
			hist->GetXaxis()->GetXbins()->GetArray(),
			hist->GetNbinsY(),
			hist->GetYaxis()->GetXbins()->GetArray());
	return cloneHist;
}

TH3D* cloneTH3D(TString name, TH3D* hist)
{
	if (!hist) return 0;
	TH3D* cloneHist = new TH3D(name,
			"",
			hist->GetNbinsX(),
			hist->GetXaxis()->GetXbins()->GetArray(),
			hist->GetNbinsY(),
			hist->GetYaxis()->GetXbins()->GetArray(),
			hist->GetNbinsZ(),
			hist->GetZaxis()->GetXbins()->GetArray());
	return cloneHist;
}

// function to make plot of 1D projection of the fit
void plotFit(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title) {
	TCanvas c1;
	c1.SetBottomMargin(0.19);
	RooPlot* plot = var.frame(min, max, nbins);
	dh->plotOn( plot );
	int iCol(0);
	int colours[6] = {kBlue, kRed, kGreen+2, kMagenta, kOrange, kGray};
	int fills[6] = {3245, 3454, 3644, 3205, 3495, 3690};

	//RooCmdArg normCmd = RooFit::NormRange("FIT");
	RooCmdArg normCmd = RooFit::Normalization(1.);

	for (std::vector<std::string>::iterator it = bkg_pdfs.begin(); it != bkg_pdfs.end(); ++it){
		if(iCol>=6) iCol=0;
		pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, RooFit::LineColor( colours[iCol] ), RooFit::LineStyle(9) );
		++iCol;
	}

	for (std::vector<std::string>::iterator it = sig_pdfs.begin(); it != sig_pdfs.end(); ++it){
		if(iCol>=6) iCol=0;
		pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, RooFit::LineColor( colours[iCol] ), RooFit::LineStyle(kDashed) );
		pdf.plotOn( plot, RooFit::Components( (*it).c_str() ), normCmd, RooFit::LineColor( colours[iCol] ), RooFit::VLines(), RooFit::FillStyle(fills[iCol]), RooFit::FillColor( colours[iCol] ), RooFit::DrawOption("LF") );
		++iCol;
	}

	pdf.plotOn( plot, normCmd, RooFit::LineColor(kBlack) );
	dh->plotOn( plot );
	plot->SetXTitle(title);
	plot->SetTitle("");
	plot->GetXaxis()->SetLabelOffset(0.02);
	plot->GetXaxis()->SetTitleOffset(1.18);
	plot->Draw();

	c1.SaveAs(savedir+"/"+name+"_fit.pdf");
	c1.SaveAs(savedir+"/"+name+"_fit.png");
	c1.SetLogy();
	c1.SaveAs(savedir+"/"+name+"_fit_log.pdf");

	delete plot;
}

// function to print parameters to a file
void printParams(TString file, const RooArgList& params) {
	FILE * pFile = fopen((savedir+"/"+file).Data(), "w");

	int nPar = params.getSize();

	for ( int i=0; i<nPar; ++i) {
		RooRealVar* par = dynamic_cast<RooRealVar*>(params.at(i));

		TString title = par->getTitle();
		float value = par->getValV();
		float error = par->getError();

		int exponent(0);
		int exponentErr(0);
		int precision(0);

		while(TMath::Abs(value) < TMath::Power(10,exponent)) exponent-=1;
		while(TMath::Abs(value) > TMath::Power(10,exponent+1)) exponent+=1;
		while(TMath::Abs(error) < TMath::Power(10,exponentErr)) exponentErr-=1;
		while(TMath::Abs(error) > TMath::Power(10,exponentErr+1)) exponentErr+=1;

		if(error < 3.5*TMath::Power(10,exponentErr)) precision = 1-exponentErr;
		else precision = -exponentErr;
		if(precision<0) precision=0;

		if(exponent<-2) {
			title+=" ($10^{";
			title+=exponent;
			title+="}$)";

			precision+=exponent;

			value/=TMath::Power(10,exponent);
			error/=TMath::Power(10,exponent);
		}

		TString format = "%-30s";
		format+="\t& $% 6.";
		format+=precision;
		format+="f \\pm % 6.";
		format+=precision;
		format+="f$ \\\\\n";

		fprintf(pFile, format, title.Data(), value, error);
	}
	fclose(pFile);
}

double getPtCorrFactor(jetType type, double ptMin, double ptMax) {
	TFile* f(0);
	TTree* t(0);
	TString cut1;
	TString cut2;
	double corr(1.);

	cut1 = "JetPT>"; cut1+=ptMin; cut1+=" && JetPT<"; cut1+=ptMax;
	cut2 = "JetPT>"; cut2+=ptMin; cut2+=" && JetPT<"; cut2+=ptMax;

	switch(type) {
		case jetRecoD04:
			std::cout << "getting pT correction factor for c->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			return corr;//no correction needed
		case jetRecoSV4:
			std::cout << "getting pT correction factor for c->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(charmSimFile);
			cut1+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && TRUEDPT>5000. && SVM[0]";
			cut2+=" && JetTRUEc && !JetTRUEb && TRUEDID[0] && SVM[0]";
			break;
		case jetRecoD05:
			std::cout << "getting pT correction factor for b->D0 in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(beautySimFile);
			//cut1+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEBPT>5000.";
			//cut2+=" && JetTRUEb && !JetTRUEc && TRUEDID[0] && TRUEBID[0] && TRUEDTRUEB>-1 && TRUEDPT>5000.";
			cut1+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEBPT[TRUEDTRUEB]>5.e3";
			cut2+=" && JetTRUEb && !JetTRUEc && TMath::Abs(TRUEDID)==421 && TRUEDTRK0IDX!=-1 && TRUEDTRK2IDX==-1 && TRUEDTRUEB!=-1 && TRUEDPT>5.e3";
			break;
		case jetRecoSV5:
			std::cout << "getting pT correction factor for b->SV in " << ptMin << " to " << ptMax << " bin..." << std::endl;
			f = TFile::Open(beautySimFile);
			cut1+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && TRUEBPT>5000. && SVM[0]";
			cut2+=" && JetTRUEb && !JetTRUEc && TRUEBID[0] && SVM[0]";
			break;
	}

	if(!f) return 1.;
	t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return 1.;

	corr = t->GetEntries(cut1)/static_cast<double>(t->GetEntries(cut2));
	std::cout << "Correction factor: " << corr << std::endl;
	f->Close();
	return corr;
}

void makeTrainTestSamples(int sample) {
	std::cout << "INFO : making samples for MC study " << sample << std::endl;
	TRandom3 rand(1000+sample);
	TString name = "_"; name+=sample;
	TFile* fin4 = TFile::Open(charmSimFile);
	TTree* tin4 = dynamic_cast<TTree*>(fin4->Get("T"));

	TFile* fout40 = TFile::Open("/tmp/dcraik/for_yandex_data_new_14X_train"+name+".root","RECREATE");
	TTree* tout40 = tin4->CloneTree(0);

	TFile* fout41 = TFile::Open("/tmp/dcraik/for_yandex_data_new_14X_test"+name+".root","RECREATE");
	TTree* tout41 = tin4->CloneTree(0);

	TFile* fin5 = TFile::Open(beautySimFile);
	TTree* tin5 = dynamic_cast<TTree*>(fin5->Get("T"));

	TFile* fout50 = TFile::Open("/tmp/dcraik/for_yandex_data_new_15X_train"+name+".root","RECREATE");
	TTree* tout50 = tin5->CloneTree(0);

	TFile* fout51 = TFile::Open("/tmp/dcraik/for_yandex_data_new_15X_test"+name+".root","RECREATE");
	TTree* tout51 = tin5->CloneTree(0);

	boost::progress_display progress( tin4->GetEntries()+tin5->GetEntries() );
	for(int ientry=0; ientry<tin4->GetEntries(); ++ientry) {
		++progress;
		tin4->GetEntry(ientry);
		if(rand.Rndm()<0.5) tout40->Fill();
		else tout41->Fill();
	}

	for(int ientry=0; ientry<tin5->GetEntries(); ++ientry) {
		++progress;
		tin5->GetEntry(ientry);
		if(rand.Rndm()<0.5) tout50->Fill();
		else tout51->Fill();
	}

	tout40->AutoSave();
	fout40->Close();
	tout41->AutoSave();
	fout41->Close();
	tout50->AutoSave();
	fout50->Close();
	tout51->AutoSave();
	fout51->Close();

	//now merge the test samples
	TChain* tout45 = new TChain("T");
	tout45->Add("/tmp/dcraik/for_yandex_data_new_14X_test"+name+".root");
	tout45->Add("/tmp/dcraik/for_yandex_data_new_15X_test"+name+".root");
	tout45->Merge("/tmp/dcraik/for_yandex_data_new"+name+".root");

	//update the standard locations and set the dataIsMC flag
	charmSimFile = "/tmp/dcraik/for_yandex_data_new_14X_train"+name+".root";
	beautySimFile = "/tmp/dcraik/for_yandex_data_new_15X_train"+name+".root";
	dataFile = "/tmp/dcraik/for_yandex_data_new"+name+".root";
	dataIsMC=true;
	dataIsResampledMC=true;
}

RooUnfoldResponse* trainUnfold(TH1* binning, jetType type, TString file) {
	TFile* f(0);
	TString nameStr;
	switch(type) {
		case jetRecoD04:
			std::cout << "INFO : training unfolding for c2D0" << std::endl;
			nameStr="_c2D0";
			break;
		case jetRecoSV4:
			std::cout << "INFO : training unfolding for c2SV" << std::endl;
			nameStr="_c2SV";
			break;
		case jetRecoD05:
			std::cout << "INFO : training unfolding for b2D0" << std::endl;
			nameStr="_b2D0";
			break;
		case jetRecoSV5:
			std::cout << "INFO : training unfolding for b2SV" << std::endl;
			nameStr="_b2SV";
			break;
		default:
			return 0;
	}

	TString fname="D0MCjets";
	if(dataIsResampledMC) fname+=file;
	fname+=nameStr+".root";

	f = TFile::Open(fname);
	if(!f) return 0;
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return 0;

	RooUnfoldResponse* resp = new RooUnfoldResponse(binning,binning);

	int countSkip(0), countMiss(0), countFake(0), countReal(0);

	double JetPT;
	double JetTruePT;
	double SVMCor;
	double D0PT;
	double weight(1.);

	t->SetBranchAddress("JetPT",       &JetPT);
	t->SetBranchAddress("JetTruePT",   &JetTruePT);

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
				if(D0PT < d0minpt ||
				   (d0maxpt!=-1 && D0PT > d0maxpt) ) continue;
				break;
			case jetRecoSV4:
				if(SVMCor<0) continue;
				break;
			case jetRecoSV5:
				if(SVMCor<0) continue;
				break;
			default:
				return 0;
		}

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

	return resp;
}

//this method doesn't try to correct acceptance and reconstruction for vertex-location effects
bool addEffsOld(TString file) {
	std::cout << "INFO : adding efficiencies to file " << file << std::endl;
	TFile* f0 = TFile::Open(dataFile);

	TFile* f2 = TFile::Open(effFile);
	TFile* f3 = TFile::Open(accFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec4 = dynamic_cast<TH2D*>(f3->Get("reff4")); //TODO this histogram has finer bins but there's currently a bug in its generation
	TH2D* hrec5 = dynamic_cast<TH2D*>(f3->Get("reff5")); //TODO this histogram has finer bins but there's currently a bug in its generation

	TH2D* hpid(0);
	if(dataIsMC) hpid = dynamic_cast<TH2D*>(f2->Get("pidmceffD045"));
	else hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	TH2D* hsel4 = dynamic_cast<TH2D*>(f2->Get("seleffD04"));
	TH2D* hsel5 = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
	//TH2D* hacc  = dynamic_cast<TH2D*>(f2->Get("acceffD045")); //superseded by RapidSim version
	//TH2D* hrec  = dynamic_cast<TH2D*>(f2->Get("receffD045"));

	if(!hacc || !hrec4 || !hrec5 || !hpid || !hsel4 || !hsel5) return false;

	std::vector<double> *vD0M = new std::vector<double>();
	std::vector<double> *vD0IPCHI2 = new std::vector<double>();
	std::vector<double> *vD0PT = new std::vector<double>();
	std::vector<double> *vD0PX = new std::vector<double>();
	std::vector<double> *vD0PY = new std::vector<double>();
	std::vector<double> *vD0PZ = new std::vector<double>();
	std::vector<double> *vD0E = new std::vector<double>();
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

	t0->SetBranchAddress("D0M",           &vD0M);
	t0->SetBranchAddress("D0IPCHI2",      &vD0IPCHI2);
	t0->SetBranchAddress("D0PT",          &vD0PT);
	t0->SetBranchAddress("D0PX",          &vD0PX);
	t0->SetBranchAddress("D0PY",          &vD0PY);
	t0->SetBranchAddress("D0PZ",          &vD0PZ);
	t0->SetBranchAddress("D0E",           &vD0E);
	t0->SetBranchAddress("D0KP",          &vD0KP);
	t0->SetBranchAddress("D0KPT",         &vD0KPT);
	t0->SetBranchAddress("D0KPX",         &vD0KPX);
	t0->SetBranchAddress("D0KPY",         &vD0KPY);
	t0->SetBranchAddress("D0KPZ",         &vD0KPZ);
	t0->SetBranchAddress("D0PIP",         &vD0PIP);
	t0->SetBranchAddress("D0PIPT",        &vD0PIPT);
	t0->SetBranchAddress("D0PIPX",        &vD0PIPX);
	t0->SetBranchAddress("D0PIPY",        &vD0PIPY);
	t0->SetBranchAddress("D0PIPZ",        &vD0PIPZ);
	t0->SetBranchAddress("D0KPNNK",       &vD0KPNNK);
	t0->SetBranchAddress("D0PIPNNPI",     &vD0PIPNNPI);

	t0->SetBranchAddress("JetPT",         &JetPT);
	t0->SetBranchAddress("JetEta",        &JetEta);
	if(dataIsMC) t0->SetBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = t0->GetEntries();

	TFile* fout = TFile::Open("D0jets"+file+".root","RECREATE");
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

	if(dataIsMC) tout->Branch("JetTruePT", &JetTruePT);

	//first make non-vector tree for fits
	boost::progress_display progress( nentries0 );
	for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
		++progress;
		t0->GetEntry(ientry);

		for(unsigned int s=0; s<vD0M->size(); ++s) {
			//first check in our tight acceptance
			TVector3 D0P (vD0PX->at(s)  ,vD0PY->at(s)  ,vD0PZ->at(s));
			TVector3 D0P0(vD0KPX->at(s) ,vD0KPY->at(s) ,vD0KPZ->at(s));
			TVector3 D0P1(vD0PIPX->at(s),vD0PIPY->at(s),vD0PIPZ->at(s));

			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>d0minpt)) continue;
			if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;

			//check PID cuts
			if(vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
			//if(vD0PIPNNPI->at(s)<0.1) continue; //pion PID turned off

			D0M         = vD0M->at(s);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
			D0PT        = D0P.Pt();
			D0Eta       = D0P.Eta();

			//get efficiency corrections
			double effacc(0.), effrec4(0.), effrec5(0.), effsel4(0.), effsel5(0.), effpid(0.);

			if(D0PT>=100000.) D0PT=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(D0PT/1000.,D0Eta));
			effrec4=      hrec4->GetBinContent(hrec4->FindBin(D0PT/1000.,D0Eta));
			effrec5=      hrec5->GetBinContent(hrec5->FindBin(D0PT/1000.,D0Eta));
			effsel4=      hsel4->GetBinContent(hsel4->FindBin(D0PT      ,D0Eta));
			effsel5=      hsel5->GetBinContent(hsel5->FindBin(D0PT      ,D0Eta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(D0PT      ,D0Eta));

			double eff4=effacc*effrec4*effsel4*effpid;
			double eff5=effacc*effrec5*effsel5*effpid;

			if(eff4<0.01 || eff5<0.01) {
//TODO//				std::cout << D0PT << "\t" << D0Eta << "\t" << effacc << "\t" << effrec4 << "\t" << effrec5 << "\t" << effsel4 << "\t" << effsel5 << "\t" << effpid << std::endl;
//TODO//				continue;
			}

			weight4 = 1./eff4;
			weight5 = 1./eff5;

			if(dataIsMC) {//scale weights for roughly continuous jet true pT
				if(JetTruePT>50000.) {
					weight4*=0.007;
					weight5*=0.007;
				} else if(JetTruePT>20000.) {
					weight4*=0.10;
					weight5*=0.10;
				} else if(JetTruePT>15000.) {
					weight4*=0.25;
					weight5*=0.25;
				}
			}

			tout->Fill();
			break;//only keep one D0 candidate per entry
		}
	}

	tout->AutoSave();
	fout->Close();

	return true;
}


bool addEffs(TString file) {
	std::cout << "INFO : adding efficiencies to file " << file << std::endl;
	TFile* f0 = TFile::Open(dataFile);

	TFile* f2 = TFile::Open(effFile);
	TFile* f3 = TFile::Open(accFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

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
	c1.SaveAs(savedir+"/D0EffAcc.pdf");
	hrec->Draw("colz");
	c1.SaveAs(savedir+"/D0EffRec.pdf");
	hcor->Draw("colz");
	c1.SaveAs(savedir+"/D0EffRecCor.pdf");
	hsel4->Draw("colz");
	c1.SaveAs(savedir+"/D0EffSelPrompt.pdf");
	hsel5->Draw("colz");
	c1.SaveAs(savedir+"/D0EffSelDispl.pdf");
	hpid->Draw("colz");
	c1.SaveAs(savedir+"/D0EffPID.pdf");

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

	t0->SetBranchAddress("D0M",           &vD0M);
	t0->SetBranchAddress("D0IPCHI2",      &vD0IPCHI2);
	t0->SetBranchAddress("D0PT",          &vD0PT);
	t0->SetBranchAddress("D0PX",          &vD0PX);
	t0->SetBranchAddress("D0PY",          &vD0PY);
	t0->SetBranchAddress("D0PZ",          &vD0PZ);
	t0->SetBranchAddress("D0E",           &vD0E);
	t0->SetBranchAddress("D0X",           &vD0X);
	t0->SetBranchAddress("D0Y",           &vD0Y);
	t0->SetBranchAddress("D0Z",           &vD0Z);
	t0->SetBranchAddress("D0KP",          &vD0KP);
	t0->SetBranchAddress("D0KPT",         &vD0KPT);
	t0->SetBranchAddress("D0KPX",         &vD0KPX);
	t0->SetBranchAddress("D0KPY",         &vD0KPY);
	t0->SetBranchAddress("D0KPZ",         &vD0KPZ);
	t0->SetBranchAddress("D0PIP",         &vD0PIP);
	t0->SetBranchAddress("D0PIPT",        &vD0PIPT);
	t0->SetBranchAddress("D0PIPX",        &vD0PIPX);
	t0->SetBranchAddress("D0PIPY",        &vD0PIPY);
	t0->SetBranchAddress("D0PIPZ",        &vD0PIPZ);
	t0->SetBranchAddress("D0KPNNK",       &vD0KPNNK);
	t0->SetBranchAddress("D0PIPNNPI",     &vD0PIPNNPI);
	t0->SetBranchAddress("D0KWEIGHT",     &vD0KWEIGHT);
	t0->SetBranchAddress("D0PIWEIGHT",    &vD0PIWEIGHT);

	t0->SetBranchAddress("JetPT",         &JetPT);
	t0->SetBranchAddress("JetEta",        &JetEta);
	if(dataIsMC) t0->SetBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = t0->GetEntries();

	TFile* fout = TFile::Open("D0jets"+file+".root","RECREATE");
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

	if(dataIsMC) tout->Branch("JetTruePT", &JetTruePT);

	//first make non-vector tree for fits
	boost::progress_display progress( nentries0 );
	for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
		++progress;
		t0->GetEntry(ientry);

		for(unsigned int s=0; s<vD0M->size(); ++s) {
			//first check in our tight acceptance
			TVector3 D0P (vD0PX->at(s)  ,vD0PY->at(s)  ,vD0PZ->at(s));
			TVector3 D0P0(vD0KPX->at(s) ,vD0KPY->at(s) ,vD0KPZ->at(s));
			TVector3 D0P1(vD0PIPX->at(s),vD0PIPY->at(s),vD0PIPZ->at(s));

			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>d0minpt)) continue;
			if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;

			//check PID cuts
			if(!dataIsMC && vD0KPNNK->at(s)<0.2 && D0P0.Pt()<25000. && D0P0.Mag()<500000.) continue; //kaon PID turned off for high p or pT
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
			if(useRhoZEffCor) effcor = hcor ->GetBinContent(hcor ->FindBin(D0RhoSq   ,D0Z  ));
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

			if(dataIsMC) {
				//put the PID efficiency into the weights for MC
				if(D0P0.Pt()<25000. && D0P0.Mag()<500000.) {
					weight4 *= vD0KWEIGHT->at(s);
					weight5 *= vD0KWEIGHT->at(s);
				}
				//currently no pion PID
				//if(D0P1.Pt()<25000. && D0P1.Mag()<500000.) {
				//	weight4 *= vD0PIWEIGHT->at(s);
				//	weight5 *= vD0PIWEIGHT->at(s);
				//}
				//scale weights for roughly continuous jet true pT
				if(JetTruePT>50000.) {
					weight4*=0.007;
					weight5*=0.007;
				} else if(JetTruePT>20000.) {
					weight4*=0.10;
					weight5*=0.10;
				} else if(JetTruePT>15000.) {
					weight4*=0.25;
					weight5*=0.25;
				}
			}

			tout->Fill();
			break;//only keep one D0 candidate per entry
		}
	}

	tout->AutoSave();
	fout->Close();

	return true;
}

bool addEffsSimple(TString file) {
	std::cout << "INFO : adding efficiencies to file " << file << std::endl;
	TFile* f0 = TFile::Open(dataFile);

	TFile* f2 = TFile::Open(simpleEffFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

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

	t0->SetBranchAddress("D0M",           &vD0M);
	t0->SetBranchAddress("D0IPCHI2",      &vD0IPCHI2);
	t0->SetBranchAddress("D0PT",          &vD0PT);
	t0->SetBranchAddress("D0PX",          &vD0PX);
	t0->SetBranchAddress("D0PY",          &vD0PY);
	t0->SetBranchAddress("D0PZ",          &vD0PZ);
	t0->SetBranchAddress("D0E",           &vD0E);
	t0->SetBranchAddress("D0X",           &vD0X);
	t0->SetBranchAddress("D0Y",           &vD0Y);
	t0->SetBranchAddress("D0Z",           &vD0Z);
	t0->SetBranchAddress("D0KP",          &vD0KP);
	t0->SetBranchAddress("D0KPT",         &vD0KPT);
	t0->SetBranchAddress("D0KPX",         &vD0KPX);
	t0->SetBranchAddress("D0KPY",         &vD0KPY);
	t0->SetBranchAddress("D0KPZ",         &vD0KPZ);
	t0->SetBranchAddress("D0PIP",         &vD0PIP);
	t0->SetBranchAddress("D0PIPT",        &vD0PIPT);
	t0->SetBranchAddress("D0PIPX",        &vD0PIPX);
	t0->SetBranchAddress("D0PIPY",        &vD0PIPY);
	t0->SetBranchAddress("D0PIPZ",        &vD0PIPZ);
	t0->SetBranchAddress("D0KPNNK",       &vD0KPNNK);
	t0->SetBranchAddress("D0PIPNNPI",     &vD0PIPNNPI);

	t0->SetBranchAddress("JetPT",         &JetPT);
	t0->SetBranchAddress("JetEta",        &JetEta);
	if(dataIsMC) t0->SetBranchAddress("JetTruePT",     &JetTruePT);

	unsigned int nentries0 = t0->GetEntries();

	TFile* fout = TFile::Open("D0jets"+file+".root","RECREATE");
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

	if(dataIsMC) tout->Branch("JetTruePT", &JetTruePT);

	//first make non-vector tree for fits
	boost::progress_display progress( nentries0 );
	for(unsigned int ientry=0; ientry<nentries0; ++ientry) {
		++progress;
		t0->GetEntry(ientry);

		for(unsigned int s=0; s<vD0M->size(); ++s) {
			//first check in our tight acceptance
			TVector3 D0P (vD0PX->at(s)  ,vD0PY->at(s)  ,vD0PZ->at(s));
			TVector3 D0P0(vD0KPX->at(s) ,vD0KPY->at(s) ,vD0KPZ->at(s));
			TVector3 D0P1(vD0PIPX->at(s),vD0PIPY->at(s),vD0PIPZ->at(s));

			if(!(D0P0.Eta()>2. && D0P0.Eta()<4.5 && D0P0.Pt()>500. && D0P0.Mag()>5000.)) continue;
			if(!(D0P1.Eta()>2. && D0P1.Eta()<4.5 && D0P1.Pt()>500. && D0P1.Mag()>5000.)) continue;
			if(!(D0P.Pt()>d0minpt)) continue;
			if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;

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

			if(dataIsMC) {//scale weights for roughly continuous jet true pT
				if(JetTruePT>50000.) {
					weight4*=0.007;
					weight5*=0.007;
				} else if(JetTruePT>20000.) {
					weight4*=0.10;
					weight5*=0.10;
				} else if(JetTruePT>15000.) {
					weight4*=0.25;
					weight5*=0.25;
				}
			}

			tout->Fill();
			break;//only keep one D0 candidate per entry
		}
	}

	tout->AutoSave();
	fout->Close();

	return true;
}

bool testEffsOld(TString file, TH1D* ptBinScheme) {
	if(!dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << file << std::endl;
	TFile* f0 = TFile::Open(dataFile);

	TFile* f2 = TFile::Open(effFile);
	TFile* f3 = TFile::Open(accFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec(0);
	if(file.BeginsWith("15")) hrec = dynamic_cast<TH2D*>(f3->Get("reff5"));
	else hrec = dynamic_cast<TH2D*>(f3->Get("reff4"));

	TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pidmceffD045"));
	TH2D* hsel(0);
	if(file.BeginsWith("15")) hsel = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
	else hsel = dynamic_cast<TH2D*>(f2->Get("seleffD04"));
	//TH2D* hacc  = dynamic_cast<TH2D*>(f2->Get("acceffD045"));
	//TH2D* hrec  = dynamic_cast<TH2D*>(f2->Get("receffD045"));

	if(!hacc || !hrec || !hpid || !hsel) return false;

	//yield histograms - 0=all; 1=passAcc; 2=passRec; 3=passSel; 4=passPID.
	//tru=truth D0; fit=eff-corrected reco D0; fnd=eff-corrected truth-matched reco D0.
	std::vector<TString> stages;
	stages.push_back("all");
	stages.push_back("geometric");
	stages.push_back("reconstruction");
	stages.push_back("selection");
	stages.push_back("PID");
	std::vector<TH1D*> truD0truePT;
	truD0truePT.push_back(cloneTH1D("tru0truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru1truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru2truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru3truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru4truept", ptBinScheme));
	std::vector<TH1D*> fitD0truePT;
	fitD0truePT.push_back(cloneTH1D("fit0truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit1truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit2truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit3truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit4truept", ptBinScheme));
	std::vector<TH1D*> fndD0truePT;
	fndD0truePT.push_back(cloneTH1D("fnd0truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd1truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd2truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd3truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd4truept", ptBinScheme));
	std::vector<TH1D*> truD0recoPT;
	truD0recoPT.push_back(cloneTH1D("tru0recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru1recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru2recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru3recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru4recopt", ptBinScheme));
	std::vector<TH1D*> fitD0recoPT;
	fitD0recoPT.push_back(cloneTH1D("fit0recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit1recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit2recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit3recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit4recopt", ptBinScheme));
	std::vector<TH1D*> fndD0recoPT;
	fndD0recoPT.push_back(cloneTH1D("fnd0recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd1recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd2recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd3recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd4recopt", ptBinScheme));

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
		fitD0truePT.at(i)->Sumw2();
		fndD0truePT.at(i)->Sumw2();
		truD0recoPT.at(i)->Sumw2();
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
					if(!(D0P.Pt()>d0minpt)) continue;
					if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;
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
			if(!(D0P.Pt()>d0minpt)) continue;
			if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;
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
			effsel =      hsel ->GetBinContent(hsel ->FindBin(dpt      ,deta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(dpt      ,deta));

			double eff=effacc*effrec*effsel*effpid;

			if(eff<0.01) {
//TODO				std::cout << dpt << "\t" << deta << "\t" << effacc << "\t" << effrec << "\t" << effsel << "\t" << effpid << std::endl;
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
			if(D0TRUEIDX->at(s)<0) break;
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
	d0Kinematics[9]->Divide(d0Kinematics[7]);
	d0Kinematics[8]->Divide(d0Kinematics[6]);
	d0Kinematics[9]->Divide(d0Kinematics[8]);
	d0Kinematics[9]->SetMinimum(0.);
	d0Kinematics[9]->SetMaximum(2.);
	d0Kinematics[9]->Draw("colz");
	c1.SaveAs(savedir+"/D0KinePIDRatio.pdf");
	d0Kinematics[7]->Divide(d0Kinematics[5]);
	d0Kinematics[6]->Divide(d0Kinematics[4]);
	d0Kinematics[7]->Divide(d0Kinematics[6]);
	d0Kinematics[7]->SetMinimum(0.);
	d0Kinematics[7]->SetMaximum(2.);
	d0Kinematics[7]->Draw("colz");
	c1.SaveAs(savedir+"/D0KineSELRatio.pdf");
	d0Kinematics[5]->Divide(d0Kinematics[3]);
	d0Kinematics[4]->Divide(d0Kinematics[2]);
	d0Kinematics[5]->Divide(d0Kinematics[4]);
	d0Kinematics[5]->SetMinimum(0.);
	d0Kinematics[5]->SetMaximum(2.);
	d0Kinematics[5]->Draw("colz");
	c1.SaveAs(savedir+"/D0KineRECRatio.pdf");
	d0Kinematics[3]->Divide(d0Kinematics[1]);
	d0Kinematics[2]->Divide(d0Kinematics[0]);
	d0Kinematics[3]->Divide(d0Kinematics[2]);
	d0Kinematics[3]->SetMinimum(0.);
	d0Kinematics[3]->SetMaximum(2.);
	d0Kinematics[3]->Draw("colz");
	c1.SaveAs(savedir+"/D0KineACCRatio.pdf");
	//for(unsigned int i=0; i<d0Kinematics.size(); ++i) {
	//	d0Kinematics[i]->Draw("colz");
	//	TString plotname("D0Kinematics");
	//	plotname+=i;
	//	plotname+=".pdf";
	//	c1.SaveAs(savedir+"/"+plotname);
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
	c1.SaveAs(savedir+"/D0KineACCPulls.pdf");
	d0KinematicsPulls[1]->Draw();
	c1.SaveAs(savedir+"/D0KineRECPulls.pdf");
	d0KinematicsPulls[2]->Draw();
	c1.SaveAs(savedir+"/D0KineSELPulls.pdf");
	d0KinematicsPulls[3]->Draw();
	c1.SaveAs(savedir+"/D0KinePIDPulls.pdf");

	return true;
}

bool testEffs(TString file, TH1D* ptBinScheme) {
	if(!dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << file << std::endl;
	TFile* f0 = TFile::Open(dataFile);

	TFile* f2 = TFile::Open(effFile);
	TFile* f3 = TFile::Open(accFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec = dynamic_cast<TH2D*>(f3->Get("reff"));
	TH2D* hcor = dynamic_cast<TH2D*>(f3->Get("corr"));

	TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	//TH2D* hpid = dynamic_cast<TH2D*>(f2->Get("pidmceffD045"));//TODO
	TH2D* hsel(0);
	if(file.BeginsWith("15")) hsel = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
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
	std::vector<TH1D*> truD0truePT;
	truD0truePT.push_back(cloneTH1D("tru0truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru1truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru2truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru3truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru4truept", ptBinScheme));
	std::vector<TH1D*> effD0truePT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effD0truePT.push_back(cloneTH1D("eff0truept", ptBinScheme));
	effD0truePT.push_back(cloneTH1D("eff1truept", ptBinScheme));
	effD0truePT.push_back(cloneTH1D("eff2truept", ptBinScheme));
	effD0truePT.push_back(cloneTH1D("eff3truept", ptBinScheme));
	effD0truePT.push_back(cloneTH1D("eff4truept", ptBinScheme));
	std::vector<TH1D*> fitD0truePT;
	fitD0truePT.push_back(cloneTH1D("fit0truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit1truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit2truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit3truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit4truept", ptBinScheme));
	std::vector<TH1D*> fndD0truePT;
	fndD0truePT.push_back(cloneTH1D("fnd0truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd1truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd2truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd3truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd4truept", ptBinScheme));
	std::vector<TH1D*> truD0recoPT;
	truD0recoPT.push_back(cloneTH1D("tru0recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru1recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru2recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru3recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru4recopt", ptBinScheme));
	std::vector<TH1D*> effD0recoPT;//each of these consists of the "true" distribution for the next stage efficiency corrected for a single stage
	effD0recoPT.push_back(cloneTH1D("eff0recopt", ptBinScheme));
	effD0recoPT.push_back(cloneTH1D("eff1recopt", ptBinScheme));
	effD0recoPT.push_back(cloneTH1D("eff2recopt", ptBinScheme));
	effD0recoPT.push_back(cloneTH1D("eff3recopt", ptBinScheme));
	effD0recoPT.push_back(cloneTH1D("eff4recopt", ptBinScheme));
	std::vector<TH1D*> fitD0recoPT;
	fitD0recoPT.push_back(cloneTH1D("fit0recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit1recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit2recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit3recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit4recopt", ptBinScheme));
	std::vector<TH1D*> fndD0recoPT;
	fndD0recoPT.push_back(cloneTH1D("fnd0recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd1recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd2recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd3recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd4recopt", ptBinScheme));

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
			if(useRhoZEffCor) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));//TODO need new datasets
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
					if(!(D0P.Pt()>d0minpt)) continue;
					if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;
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
			if(!(D0P.Pt()>d0minpt)) continue;
			if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;
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
			if(useRhoZEffCor) effrec*= hcor ->GetBinContent(hcor ->FindBin(rhoSq    ,z));
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
	c1.SaveAs(savedir+"/D0KinePIDRatio.pdf");
	d0Kinematics[7]->Divide(d0Kinematics[5]);
	d0Kinematics[6]->Divide(d0Kinematics[4]);
	d0Kinematics[7]->Divide(d0Kinematics[6]);
	d0Kinematics[7]->SetMinimum(0.);
	d0Kinematics[7]->SetMaximum(2.);
	d0Kinematics[7]->Draw("colz");
	c1.SaveAs(savedir+"/D0KineSELRatio.pdf");
	d0Kinematics[5]->Divide(d0Kinematics[3]);
	d0Kinematics[4]->Divide(d0Kinematics[2]);
	d0Kinematics[5]->Divide(d0Kinematics[4]);
	d0Kinematics[5]->SetMinimum(0.);
	d0Kinematics[5]->SetMaximum(2.);
	d0Kinematics[5]->Draw("colz");
	c1.SaveAs(savedir+"/D0KineRECRatio.pdf");
	d0Kinematics[3]->Divide(d0Kinematics[1]);
	d0Kinematics[2]->Divide(d0Kinematics[0]);
	d0Kinematics[3]->Divide(d0Kinematics[2]);
	d0Kinematics[3]->SetMinimum(0.);
	d0Kinematics[3]->SetMaximum(2.);
	d0Kinematics[3]->Draw("colz");
	c1.SaveAs(savedir+"/D0KineACCRatio.pdf");
	//for(unsigned int i=0; i<d0Kinematics.size(); ++i) {
	//	d0Kinematics[i]->Draw("colz");
	//	TString plotname("D0Kinematics");
	//	plotname+=i;
	//	plotname+=".pdf";
	//	c1.SaveAs(savedir+"/"+plotname);
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
	c1.SaveAs(savedir+"/D0KineACCPulls.pdf");
	d0KinematicsPulls[1]->Draw();
	c1.SaveAs(savedir+"/D0KineRECPulls.pdf");
	d0KinematicsPulls[2]->Draw();
	c1.SaveAs(savedir+"/D0KineSELPulls.pdf");
	d0KinematicsPulls[3]->Draw();
	c1.SaveAs(savedir+"/D0KinePIDPulls.pdf");

	return true;
}

bool testEffsSimple(TString file, TH1D* ptBinScheme) {
	if(!dataIsMC) return false;
	std::cout << "INFO : testing efficiencies for file " << file << std::endl;
	TFile* f0 = TFile::Open(dataFile);

	TFile* f2 = TFile::Open(simpleEffFile);

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* heff(0);
	if(file.BeginsWith("15")) heff = dynamic_cast<TH2D*>(f2->Get("effD05"));
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
	truD0truePT.push_back(cloneTH1D("tru0truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru1truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru2truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru3truept", ptBinScheme));
	truD0truePT.push_back(cloneTH1D("tru4truept", ptBinScheme));
	std::vector<TH1D*> fitD0truePT;
	fitD0truePT.push_back(cloneTH1D("fit0truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit1truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit2truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit3truept", ptBinScheme));
	fitD0truePT.push_back(cloneTH1D("fit4truept", ptBinScheme));
	std::vector<TH1D*> fndD0truePT;
	fndD0truePT.push_back(cloneTH1D("fnd0truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd1truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd2truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd3truept", ptBinScheme));
	fndD0truePT.push_back(cloneTH1D("fnd4truept", ptBinScheme));
	std::vector<TH1D*> truD0recoPT;
	truD0recoPT.push_back(cloneTH1D("tru0recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru1recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru2recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru3recopt", ptBinScheme));
	truD0recoPT.push_back(cloneTH1D("tru4recopt", ptBinScheme));
	std::vector<TH1D*> fitD0recoPT;
	fitD0recoPT.push_back(cloneTH1D("fit0recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit1recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit2recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit3recopt", ptBinScheme));
	fitD0recoPT.push_back(cloneTH1D("fit4recopt", ptBinScheme));
	std::vector<TH1D*> fndD0recoPT;
	fndD0recoPT.push_back(cloneTH1D("fnd0recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd1recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd2recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd3recopt", ptBinScheme));
	fndD0recoPT.push_back(cloneTH1D("fnd4recopt", ptBinScheme));

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
	for(int i=0; i< ptBinScheme->GetNbinsX(); ++i) {
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

		int ptBin = ptBinScheme->FindBin(JetPT);
		if(ptBin<1||ptBin>ptBinScheme->GetNbinsX()) continue;

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
					if(!(D0P.Pt()>d0minpt)) continue;
					if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;
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
			if(!(D0P.Pt()>d0minpt)) continue;
			if(!(d0maxpt==-1 || D0P.Pt() < d0maxpt) ) continue;
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
	c1.SaveAs(savedir+"/D0KineRatio.pdf");
	c1.SetLogx(0);
	
	for(int i=0; i< ptBinScheme->GetNbinsX(); ++i) {
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
		c1.SaveAs(savedir+"/D0Mass"+suffix+".pdf");
	}
	
	return true;
}

bool weightMC(TString file, jetType type) {
	TString nameStr="";

	TFile* fin(0);
	switch(type) {
		case jetRecoD04:
			std::cout << "INFO : weighting c2D0 MC sample for job " << file << std::endl;
			nameStr="_c2D0";
			fin = TFile::Open(charmSimFile);
			break;
		case jetRecoSV4:
			std::cout << "INFO : weighting c2SV MC sample for job " << file << std::endl;
			nameStr="_c2SV";
			fin = TFile::Open(charmSimFile);
			break;
		case jetRecoD05:
			std::cout << "INFO : weighting b2D0 MC sample for job " << file << std::endl;
			nameStr="_b2D0";
			fin = TFile::Open(beautySimFile);
			break;
		case jetRecoSV5:
			std::cout << "INFO : weighting b2SV MC sample for job " << file << std::endl;
			nameStr="_b2SV";
			fin = TFile::Open(beautySimFile);
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

	tin->SetBranchAddress("SVMCor",      &vSVMCor);
	tin->SetBranchAddress("SVN",         &vSVN);

	tin->SetBranchAddress("D0M",         &vD0M);
	tin->SetBranchAddress("D0PT",        &vD0PT);
	tin->SetBranchAddress("D0IPCHI2",    &vD0IPCHI2);
	tin->SetBranchAddress("D0KWEIGHT",   &vD0KWEIGHT);
	tin->SetBranchAddress("D0PIWEIGHT",  &vD0PIWEIGHT);

	//weight MC for continuous true jet pT
	unsigned int npt=4;
	double* ptInputWeightBins  = new double[npt +1]{10000.,15000.,20000.,50000.,100000.};

	TH1D* jetTruePtWeights = new TH1D("jetTruePtWeights","",npt,ptInputWeightBins);
	if(type==jetRecoD04 || type==jetRecoSV4) {
		jetTruePtWeights->SetBinContent(1,1.);
		jetTruePtWeights->SetBinContent(2,0.3);
		jetTruePtWeights->SetBinContent(3,0.12);
		jetTruePtWeights->SetBinContent(4,0.008);
	} else if(type==jetRecoD05 || type==jetRecoSV5) {
		jetTruePtWeights->SetBinContent(1,1.);
		jetTruePtWeights->SetBinContent(2,0.2);
		jetTruePtWeights->SetBinContent(3,0.08);
		jetTruePtWeights->SetBinContent(4,0.006);
	} else {
		return false;
	}
	//end true(pT) weights setup

	//for D0 also weight to match pT(D0)/pT(jet) to data
	TFile* fdata = TFile::Open("D0jets"+file+".root");
	TTree* tdata = dynamic_cast<TTree*>(fdata->Get("T"));

	TH1D* fPtDWeights = new TH1D("fPtDWeights","",20,0.,2.0);

	if(type==jetRecoD04 || type==jetRecoD05) {
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
		c.SaveAs(savedir+"/D0fPtReweighting"+file+nameStr+".pdf");

		fPtDWeights->Draw();
		c.SaveAs(savedir+"/D0fPtWeights"+file+nameStr+".pdf");
	}
	//end f(pT) weights setup

	//now make the weighted output file
	TString fname="D0MCjets";
	if(dataIsResampledMC) fname+=file;
	fname+=nameStr+".root";

	TFile* fout = TFile::Open(fname,"RECREATE");
	TTree* tout = new TTree("T","");
	tout->SetDirectory(fout);

	double D0M(0.), D0PT(0.), D0LogIPChi2(0.), SVMCor(0.), SVN(0.), weight(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetTruePT",    &JetTruePT);
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

		//only use matching jets
		switch(type) {
			case jetRecoD04:
			case jetRecoD05:
				if(!JetTrueD0) continue;
				if(D0PT < d0minpt || (d0maxpt!=-1 && D0PT > d0maxpt) ) continue;
				break;
			case jetRecoSV4:
				if(!JetTrueDSV) continue;
				if(SVMCor<0) continue;
				break;
			case jetRecoSV5:
				if(!JetTrueBSV) continue;
				if(SVMCor<0) continue;
				break;
			default:
				return 0;
		}

		weight = jetTruePtWeights->GetBinContent(jetTruePtWeights->FindBin(JetTruePT));

		if(type==jetRecoD04||type==jetRecoD05) {
			weight *= fPtDWeights->GetBinContent(fPtDWeights->FindBin(D0PT/JetPT));

			//apply PID weights
			weight *= vD0KWEIGHT->at(0);
			weight *= vD0PIWEIGHT->at(0);
		}

		tout->Fill();
	}

	tout->SetDirectory(fout);
	tout->Write();
	tout->AutoSave();
	fout->Close();

	return true;
}

bool fitD0_1D1D(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
	std::cout << "INFO : fitting D0 - 1D mass and IP fits to file " << file << ", pT in (" << minpt << "," << maxpt << ")" << std::endl;
	bool fix(true);

	double valPromptMean(0.), valPromptWidth(0.), valPromptAsym(0.), valPromptRhoL(0.), valPromptRhoR(0.),  valDisplacedMean(5.), valDisplacedWidth(1.5);
	//D0 settings
	valPromptMean = 0.87;
	valPromptWidth = 1.09;
	valPromptAsym = -0.29;
	valPromptRhoL = 1.30;
	valPromptRhoR = 1.69;

	TString ptStr;
	ptStr+=minpt; ptStr+="-"; ptStr+=maxpt;

	if(d0minpt!=5000. || d0maxpt!=-1) {
		ptStr+="_D0";
		ptStr+=d0minpt;
		ptStr+="-";
		ptStr+=d0maxpt;
	}

	TFile* f0 = TFile::Open("D0jets"+file+".root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TString weightName;

	if(flavour==4) {
		weightName = "weight4";
	} else if(flavour==5) {
		weightName = "weight5";
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
	RooRealVar JetPT( "JetPT",  "JetPT",  minpt,  maxpt,  "MeV/#it{c}"); 
	RooRealVar DWeight(  weightName,  weightName,  0.,   1000.); 

	RooArgSet obs;
	obs.add(JetPT);
	obs.add(DM);
	obs.add(DPT);
	obs.add(DLOGIPCHI2);
	obs.add(DWeight);

	RooDataSet* ds = new RooDataSet("ds","ds", obs, RooFit::WeightVar(weightName), RooFit::Import(*t0));

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
	gSystem->RedirectOutput(savedir+"/"+dType+"_"+file+"_"+ptStr+"_fits.log","w");
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

	gSystem->RedirectOutput(savedir+"/"+dType+"_"+file+"_"+ptStr+"_fits.log","a");
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
	if(flavour==5) {
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

	plotFit(DM, massVal-massScaleFact*80., massVal+massScaleFact*80., 40, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_M_"+ptStr, "m");
	plotFit(DLOGIPCHI2, -5., 15., 20, dsPeak, data_pdf2, sig_pdfs2, bkg_pdfs2, dType+"_"+file+"_IPChi2_"+ptStr, "log(IP#chi^{2})");

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
	printParams(dType+"_"+file+"_params_"+ptStr+".dat",params);

	return true;
}

bool fitD0_SBS(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
	std::cout << "INFO : fitting D0 - sideband subtraction and IP cut to file " << file << ", pT in (" << minpt << "," << maxpt << ")"  << std::endl;
	TString ptStr;
	ptStr+=minpt; ptStr+="-"; ptStr+=maxpt;

	if(d0minpt!=5000. || d0maxpt!=-1) {
		ptStr+="_D0";
		ptStr+=d0minpt;
		ptStr+="-";
		ptStr+=d0maxpt;
	}

	TFile* f0 = TFile::Open("D0jets"+file+".root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TString weightName;

	if(flavour==4) {
		weightName = "weight4";
	} else if(flavour==5) {
		weightName = "weight5";
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
	jetPTCut+="JetPT>"; jetPTCut+=minpt; jetPTCut+=" && ";
	jetPTCut+="JetPT<"; jetPTCut+=maxpt; jetPTCut+=")";

	TString peakCut=weightName;
	peakCut+=" * (TMath::Abs(";
	peakCut+=dType;
	peakCut+="M-";
	peakCut+=massVal;
	peakCut+=")<";
	peakCut+=massScaleFact*30.;
	peakCut+=") * ";
	peakCut+=jetPTCut;
	TString sideCut=weightName;
	sideCut+=" * (TMath::Abs(";
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
	c.SaveAs(savedir+"/"+dType+"_"+file+"_IPChi2_SBS_"+ptStr+".pdf");

	std::cout << hsub.Integral() << "\t" << hsub.Integral(1,14) << "\t" << hsub.Integral(15,40) << std::endl;

	//extract results - not yet corrected for prompt<->displaced migration
	yield=0.;
	error=0.;
	if(flavour==5) {
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

bool fitD0_SBS1D(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
	std::cout << "INFO : fitting D0 - sideband subtraction and IP fit to file " << file << ", pT in (" << minpt << "," << maxpt << ")"  << std::endl;
	TString ptStr;
	ptStr+=minpt; ptStr+="-"; ptStr+=maxpt;

	if(d0minpt!=5000. || d0maxpt!=-1) {
		ptStr+="_D0";
		ptStr+=d0minpt;
		ptStr+="-";
		ptStr+=d0maxpt;
	}

	TFile* f0 = TFile::Open("D0jets"+file+".root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TString weightName;

	if(flavour==4) {
		weightName = "weight4";
	} else if(flavour==5) {
		weightName = "weight5";
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
	jetPTCut+="JetPT>"; jetPTCut+=minpt; jetPTCut+=" && ";
	jetPTCut+="JetPT<"; jetPTCut+=maxpt; jetPTCut+=")";

	TString peakCut=weightName;
	peakCut+=" * (TMath::Abs(";
	peakCut+=dType;
	peakCut+="M-";
	peakCut+=massVal;
	peakCut+=")<";
	peakCut+=massScaleFact*40.;
	peakCut+=") * ";
	peakCut+=jetPTCut;
	TString sideCut=weightName;
	sideCut+=" * (TMath::Abs(";
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

	t0->Draw(dType+"LogIPChi2>>hpeak",peakCut);
	t0->Draw(dType+"LogIPChi2>>hside",sideCut);

	t0->Draw(dType+"M>>mpeak",peakCut);
	t0->Draw(dType+"M>>mside",sideCut);
	t0->Draw(dType+"M>>mall");

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
	mall.Draw("same");
	c.SaveAs(savedir+"/"+dType+"_"+file+"_M_SBS_"+ptStr+".pdf");

	hpeak.GetXaxis()->SetTitle("log #it{#chi}^{2}_{IP}(#it{D}^{0})");
	hpeak.GetYaxis()->SetTitle("Events / (0.5)");
	hpeak.SetLineColor(kBlue);
	hpeak.Draw("E1");
	hsub.SetLineColor(kGreen+2);
	hsub.Draw("E1 same");
	hside.SetLineColor(kRed);
	hside.Draw("E1 same");
	c.SaveAs(savedir+"/"+dType+"_"+file+"_IPChi2_SBS_"+ptStr+".pdf");

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

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(prompt,displaced), RooArgList(promptYield,displacedYield) );

	gSystem->RedirectOutput(savedir+"/"+dType+"_"+file+"_"+ptStr+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));
	gSystem->RedirectOutput(0);

	//extract results - corrected for the mass window cut
	if(flavour==5) {
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

	plotFit(DLOGIPCHI2, -5., 15., 20, &dh, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_IPChi2_SBS1D_"+ptStr, "log(IP#chi^{2})");

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
	printParams(dType+"_"+file+"_params_SBS1D_"+ptStr+".dat",params);

	return true;
}

bool fitD0_2D(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
	std::cout << "INFO : fitting D0 - 2D mass and IP fit to file " << file << ", pT in (" << minpt << "," << maxpt << ")"  << std::endl;
	bool fix(true);

	double valPromptMean(0.), valPromptWidth(0.), valPromptAsym(0.), valPromptRhoL(0.), valPromptRhoR(0.),  valDisplacedMean(5.), valDisplacedWidth(1.5);
	//D0 settings
	valPromptMean = 0.87;
	valPromptWidth = 1.09;
	valPromptAsym = -0.29;
	valPromptRhoL = 1.30;
	valPromptRhoR = 1.69;

	TString ptStr;
	ptStr+=minpt; ptStr+="-"; ptStr+=maxpt;

	if(d0minpt!=5000. || d0maxpt!=-1) {
		ptStr+="_D0";
		ptStr+=d0minpt;
		ptStr+="-";
		ptStr+=d0maxpt;
	}

	TFile* f0 = TFile::Open("D0jets"+file+".root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f1 = TFile::Open("../skim/combHist.root");
	TH2D* combhist = dynamic_cast<TH2D*>(f1->Get("comb"));

	if(!t0 || !combhist) return false;

	TString weightName;

	if(flavour==4) {
		weightName = "weight4";
	} else if(flavour==5) {
		weightName = "weight5";
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
	RooRealVar JetPT( "JetPT",  "JetPT",  minpt,  maxpt,  "MeV/#it{c}"); 
	RooRealVar DWeight(  weightName,  weightName,  0.,   1000.); 

	RooArgSet obs;
	obs.add(JetPT);
	obs.add(DM);
	obs.add(DPT);
	obs.add(DLOGIPCHI2);
	obs.add(DWeight);

	RooDataSet* ds = new RooDataSet("ds","ds", obs, RooFit::WeightVar(weightName), RooFit::Import(*t0));

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

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(prompt,displaced,bkg), RooArgList(promptYield,displacedYield,bkgYield) );

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(savedir+"/"+dType+"_"+file+"_"+ptStr+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));
	gSystem->RedirectOutput(0);

	DM.setRange("signal",massVal-massScaleFact*20.,massVal+massScaleFact*20.);
	DM.setRange("sideLo",massVal-massScaleFact*80.,massVal-massScaleFact*60.);
	DM.setRange("sideHi",massVal+massScaleFact*60.,massVal+massScaleFact*80.);
	DM.setRange("full",  massVal-massScaleFact*80.,massVal+massScaleFact*80.);

	DLOGIPCHI2.setRange("sigIPCHI2",-5.,3.);
	DLOGIPCHI2.setRange("fullIPCHI2",-5.,15.);

	//extract results - corrected for the mass window cut
	if(flavour==5) {
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

	plotFit(DM, massVal-massScaleFact*80., massVal+massScaleFact*80., 40, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_M_2D_"+ptStr, "m");
	plotFit(DLOGIPCHI2, -5., 15., 20, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_IPChi2_2D_"+ptStr, "log(IP#chi^{2})");

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
	printParams(dType+"_"+file+"_params2D_"+ptStr+".dat",params);

	return true;
}

bool fitD0(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
	switch(whichFit) {
		case fit1D1D:
			return fitD0_1D1D(flavour, minpt, maxpt, yield, error, file);
		case fit2D:
			return fitD0_2D(flavour, minpt, maxpt, yield, error, file);
		case fitSidebandSub:
			return fitD0_SBS(flavour, minpt, maxpt, yield, error, file);
		case fitSBS1D:
			return fitD0_SBS1D(flavour, minpt, maxpt, yield, error, file);
		default:
			return false;
	}
}

bool getD0Yields(TH1D* hist4, TH1D* hist5, TString file) {
	double yield(0.), error(0.);

	if(!hist4 || !hist5) return false;

	for(int i=1; i<=hist4->GetNbinsX(); ++i) {
		if(!fitD0(4,hist4->GetBinLowEdge(i),hist4->GetBinLowEdge(i+1),yield,error,file)) return false;
		double corr = getPtCorrFactor(jetRecoD04,hist4->GetBinLowEdge(i),hist4->GetBinLowEdge(i+1));
		hist4->SetBinContent(i,yield*corr);
		hist4->SetBinError(  i,error*corr);
	}

	for(int i=1; i<=hist5->GetNbinsX(); ++i) {
		if(!fitD0(5,hist5->GetBinLowEdge(i),hist5->GetBinLowEdge(i+1),yield,error,file)) return false;
		double corr = getPtCorrFactor(jetRecoD05,hist5->GetBinLowEdge(i),hist5->GetBinLowEdge(i+1));
		hist5->SetBinContent(i,yield*corr);
		hist5->SetBinError(  i,error*corr);
	}

	return true;
}

//function to fit features for a single sample
TString fitSV(double& NC, double& eC, double& NB, double& eB, TString sample="100", double minPT=20000, double maxPT=30000){
	std::cout << "INFO : fitting SV - 1D corrected mass and Ntrk fit to file " << sample << ", pT in (" << minPT << "," << maxPT << ")" << std::endl;

	TString ptStr;
	ptStr+=minPT; ptStr+="-"; ptStr+=maxPT;

	TString ptCutStr;
	ptCutStr = "JetPT>"; ptCutStr+=minPT; ptCutStr+=" && JetPT<"; ptCutStr+=maxPT;

	TFile* f0 = TFile::Open(lightSimFile);
	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f4 = TFile::Open(charmSimFile);
	TTree* t4 = dynamic_cast<TTree*>(f4->Get("T"));

	TFile* f5 = TFile::Open(beautySimFile);
	TTree* t5 = dynamic_cast<TTree*>(f5->Get("T"));

	TFile* fd = TFile::Open(dataFile);
	TTree* td = dynamic_cast<TTree*>(fd->Get("T"));


	if(!t0 || !t4 || !t5 || !td) return "FAIL";

	double mmin(500.), mmax(10000.);
	int nmbins(96.);
	//int nmbins(20.);
	int nbin(nmbins+4);
	double scaleNtrk=4./nmbins;//scale down Ntrks part for nicer plots
	double min(0.), max(nbin);
	double scale(1./(1.+scaleNtrk));//scale to correct for ntrks part when extracting yields
	if(whichSVFit==fitMCorr) {
		nbin=nmbins;
		max=nmbins;
		scale=1.;
	} else if(whichSVFit==fitNTrk) {
		nbin=4;
		min=nmbins;
		scale=1.;
		scaleNtrk=1.;
	}

	TH1D svmcorsvn_0("svmcorsvn_0","",nbin,min,max);
	TH1D svmcorsvn_4("svmcorsvn_4","",nbin,min,max);
	TH1D svmcorsvn_5("svmcorsvn_5","",nbin,min,max);
	TH1D svmcorsvn_d("svmcorsvn_d","",nbin,min,max);

	svmcorsvn_4.Sumw2();
	svmcorsvn_5.Sumw2();
	svmcorsvn_d.Sumw2();

	std::vector<double>* SVN = new std::vector<double>();
	std::vector<double>* SVMCor = new std::vector<double>();
	double JetPT, JetTruePT;

	t0->SetBranchAddress("SVMCor", &SVMCor);
	t4->SetBranchAddress("SVMCor", &SVMCor);
	t5->SetBranchAddress("SVMCor", &SVMCor);
	td->SetBranchAddress("SVMCor", &SVMCor);
	t0->SetBranchAddress("SVN", &SVN);
	t4->SetBranchAddress("SVN", &SVN);
	t5->SetBranchAddress("SVN", &SVN);
	td->SetBranchAddress("SVN", &SVN);
	t0->SetBranchAddress("JetPT", &JetPT);
	t4->SetBranchAddress("JetPT", &JetPT);
	t5->SetBranchAddress("JetPT", &JetPT);
	td->SetBranchAddress("JetPT", &JetPT);
	t4->SetBranchAddress("JetTruePT", &JetTruePT);
	t5->SetBranchAddress("JetTruePT", &JetTruePT);
	if(dataIsMC) td->SetBranchAddress("JetTruePT", &JetTruePT);

	//get MC/backwards efficiencies for estimating total yields
	int eff0denom(0.), eff4denom(0.), eff5denom(0.);
	int eff0num(0.), eff4num(0.), eff5num(0.);
	double weight(1.);

	unsigned int npt=4;
	double* ptInputWeightBins  = new double[npt +1]{10000.,15000.,20000.,50000.,100000.};

	TH1D* mcWeights4 = new TH1D("mcWeights4","",npt,ptInputWeightBins);
	TH1D* mcWeights5 = new TH1D("mcWeights5","",npt,ptInputWeightBins);
	TH1D* mcWeightsD = new TH1D("mcWeightsD","",npt,ptInputWeightBins);//average to use for mixed sample
	mcWeights4->SetBinContent(1, 1.);
	mcWeights4->SetBinContent(2, 0.3);
	mcWeights4->SetBinContent(3, 0.12);
	mcWeights4->SetBinContent(4, 0.008);
	mcWeights5->SetBinContent(1, 1.);
	mcWeights5->SetBinContent(2, 0.2);
	mcWeights5->SetBinContent(3, 0.08);
	mcWeights5->SetBinContent(4, 0.006);
	mcWeightsD->SetBinContent(1, 1.);
	mcWeightsD->SetBinContent(2, 0.25);
	mcWeightsD->SetBinContent(3, 0.10);
	mcWeightsD->SetBinContent(4, 0.007);

	for(int i=0; i<t0->GetEntries(); ++i) {
		t0->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff0denom;
		if(SVN->size()<1) continue;
		++eff0num;
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_0.Fill(SVN->at(0)+nmbins-2,scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_0.Fill(nmbins+3,scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_0.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin));
		else if(SVMCor->at(0)>mmax) svmcorsvn_0.Fill(nmbins-1);
	}
	std::cout << svmcorsvn_0.Integral(1,nmbins) << "\t" << svmcorsvn_0.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	for(int i=0; i<t4->GetEntries(); ++i) {
		t4->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff4denom;
		if(SVN->size()<1) continue;
		++eff4num;
		if(JetTruePT>=100000.) JetTruePT=99999.;
		weight = mcWeights4->GetBinContent(mcWeights4->FindBin(JetTruePT));
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_4.Fill(SVN->at(0)+nmbins-2,weight*scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_4.Fill(nmbins+3,weight*scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_4.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin),weight);
		else if(SVMCor->at(0)>mmax) svmcorsvn_4.Fill(nmbins-1,weight);
	}
	std::cout << svmcorsvn_4.Integral(1,nmbins) << "\t" << svmcorsvn_4.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	for(int i=0; i<t5->GetEntries(); ++i) {
		t5->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff5denom;
		if(SVN->size()<1) continue;
		++eff5num;
		if(JetTruePT>=100000.) JetTruePT=99999.;
		weight = mcWeights5->GetBinContent(mcWeights5->FindBin(JetTruePT));
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_5.Fill(SVN->at(0)+nmbins-2,weight*scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_5.Fill(nmbins+3,weight*scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_5.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin),weight);
		else if(SVMCor->at(0)>mmax) svmcorsvn_5.Fill(nmbins-1,weight);
	}
	std::cout << svmcorsvn_5.Integral(1,nmbins) << "\t" << svmcorsvn_5.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	for(int i=0; i<td->GetEntries(); ++i) {
		weight=1.;
		if(dataIsMC) weight = mcWeightsD->GetBinContent(mcWeightsD->FindBin(JetTruePT));
		td->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		if(SVN->size()<1) continue;
		if(SVN->at(0)>1.5 && SVN->at(0)<4.5) svmcorsvn_d.Fill(SVN->at(0)+nmbins-2,weight*scaleNtrk);
		else if(SVN->at(0)>4.5) svmcorsvn_d.Fill(nmbins+3,weight*scaleNtrk);
		if(SVMCor->at(0)>mmin && SVMCor->at(0)<mmax) svmcorsvn_d.Fill((nmbins-1)*(SVMCor->at(0)-mmin)/(mmax-mmin),weight);
		else if(SVMCor->at(0)>mmax) svmcorsvn_d.Fill(nmbins-1,weight);
	}
	std::cout << svmcorsvn_d.Integral(1,nmbins) << "\t" << svmcorsvn_d.Integral(nmbins+1,nbin)/scaleNtrk << std::endl;

	// -- variables from datasets
	RooRealVar SVComb(  "SVComb",  "SVComb",  min, max,  ""); 
	SVComb.setBins(nbin);

	// -- variables for shape
	RooRealVar yieldB(  "yieldB",   "yieldB",  0.4*svmcorsvn_d.Integral(), 0., svmcorsvn_d.Integral());
	RooRealVar yieldC(  "yieldC",   "yieldC",  0.3*svmcorsvn_d.Integral(), 0., svmcorsvn_d.Integral());
	RooRealVar yieldQ(  "yieldQ",   "yieldQ",  0.3*svmcorsvn_d.Integral(), 0., svmcorsvn_d.Integral());

	RooDataHist histSVMSVNB("histCombB", "histCombB", RooArgList(SVComb), &svmcorsvn_5);
	RooDataHist histSVMSVNC("histCombC", "histCombC", RooArgList(SVComb), &svmcorsvn_4);
	RooDataHist histSVMSVNQ("histCombQ", "histCombQ", RooArgList(SVComb), &svmcorsvn_0);

	// -- simulation PDFs for each category
	RooHistPdf pdfB( "pdfB", "pdfB", RooArgSet(SVComb), histSVMSVNB ); 
	RooHistPdf pdfC( "pdfC", "pdfC", RooArgSet(SVComb), histSVMSVNC ); 
	RooHistPdf pdfQ( "pdfQ", "pdfQ", RooArgSet(SVComb), histSVMSVNQ ); 

	// -- total PDFs (with and without normalisation)
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(pdfB, pdfC, pdfQ), RooArgList(yieldB, yieldC, yieldQ) );

	// -- add all feature observables to dataset
	RooArgSet obs;
	obs.add(SVComb);

	RooDataHist dh("dh", "dh", RooArgList(SVComb), &svmcorsvn_d);

	// -- fit model pdf to the dataset ----------------------------------------------
	gSystem->RedirectOutput(savedir+"/SVComb_"+sample+"_"+ptStr+"_fits.log","w");
	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));
	gSystem->RedirectOutput(0);
	////TODO fudge c+10% from b
	//double fudgeShift = yieldC.getValV()*0.1;
	//yieldC.setVal(yieldC.getVal()+fudgeShift);
	//yieldB.setVal(yieldB.getVal()-fudgeShift);

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ" );
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	std::vector<std::string> bkg_pdfs;

	TString plotName = "SVComb_"+sample+"_"; plotName+=ptStr;
	plotFit(SVComb, min, max, nbin, &dh, data_pdf, sig_pdfs, bkg_pdfs, plotName, "M_{cor}, N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_comb_"+sample+"_"; paramsName+=minPT; paramsName+="-"; paramsName+=maxPT; paramsName+=".dat";
	printParams(paramsName,params);

	double NQ, Ntot;
	double eQ;

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV()*scale;
	NC = yieldC.getValV()*scale;
	NQ = yieldQ.getValV()*scale;
	eB = yieldB.getError()*scale;
	eC = yieldC.getError()*scale;
	eQ = yieldQ.getError()*scale;
	Ntot*=scale;
	double eff5 = static_cast<double>(eff5num)/eff5denom;
	double eff4 = static_cast<double>(eff4num)/eff4denom;
	double eff0 = static_cast<double>(eff0num)/eff0denom;
	double erreff5 = TMath::Sqrt(static_cast<double>(eff5num))/eff5denom;
	double erreff4 = TMath::Sqrt(static_cast<double>(eff4num))/eff4denom;
	double erreff0 = TMath::Sqrt(static_cast<double>(eff0num))/eff0denom;
	double corr5 = NB/eff5;
	double corr4 = NC/eff4;
	double corr0 = NQ/eff0;
	double errcorr5 = corr5 * TMath::Sqrt(eB*eB/NB/NB + erreff5*erreff5/eff5/eff5);
	double errcorr4 = corr4 * TMath::Sqrt(eC*eC/NC/NC + erreff4*erreff4/eff4/eff4);
	double errcorr0 = corr0 * TMath::Sqrt(eQ*eQ/NQ/NQ + erreff0*erreff0/eff0/eff0);
	printf("% 6.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\t% 6.0f+/-% 4.0f\n\n", Ntot, NB, eB, NC, eC, NQ, eQ, corr5, errcorr5, corr4, errcorr4, corr0, errcorr0);

	return true;
}

bool getSVYields(TH1D* hist4, TH1D* hist5, TString file) {
	double nB(0.), eB(0.), nC(0.), eC(0.);

	if(!hist4 || !hist5) return false;

	for(int i=1; i<=hist4->GetNbinsX(); ++i) {
		if(!fitSV(nC,eC,nB,eB,file,hist4->GetBinLowEdge(i),hist4->GetBinLowEdge(i+1))) return false;
		double corrC = getPtCorrFactor(jetRecoSV4,hist4->GetBinLowEdge(i),hist4->GetBinLowEdge(i+1));
		double corrB = getPtCorrFactor(jetRecoSV5,hist5->GetBinLowEdge(i),hist5->GetBinLowEdge(i+1));
		hist4->SetBinContent(i,nC*corrC);
		hist4->SetBinError(  i,eC*corrC);
		hist5->SetBinContent(i,nB*corrB);
		hist5->SetBinError(  i,eB*corrB);
	}

	return true;
}

bool getTruth(TH1D* trueD04, TH1D* trueD05, TH1D* trueSV4, TH1D* trueSV5, bool useTruePT=false) {
	std::cout << "INFO : getting truth information" << std::endl;
	if(!trueD04 || !trueD05 || !trueSV4 || !trueSV5) return false;

	TFile* f = TFile::Open(dataFile);
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

	double JetPT;
	double JetTruePT;
	double JetTrueD0;
	double JetTrueDSV;
	double JetTrueBSV;

	std::vector<double>* TRUEDID = new std::vector<double>();
	std::vector<double>* TRUEDFROMB = new std::vector<double>();

	std::vector<double>* SVN = new std::vector<double>();

	t->SetBranchAddress("JetPT",      &JetPT);
	t->SetBranchAddress("JetTruePT",  &JetTruePT);
	t->SetBranchAddress("JetTRUED0",  &JetTrueD0);
	t->SetBranchAddress("JetTRUEDSV", &JetTrueDSV);
	t->SetBranchAddress("JetTRUEBSV", &JetTrueBSV);

	t->SetBranchAddress("TRUEDID",    &TRUEDID);
	t->SetBranchAddress("TRUEDFROMB", &TRUEDFROMB);

	t->SetBranchAddress("SVN",        &SVN);

	boost::progress_display progress( t->GetEntries() );
	for(int i=0; i<t->GetEntries(); ++i) {
		++progress;
		t->GetEntry(i);
		if(useTruePT) JetPT = JetTruePT;

		double weight(1.);
		if(JetTruePT>50000.) {
			weight=0.007;
		} else if(JetTruePT>20000.) {
			weight=0.10;
		} else if(JetTruePT>15000.) {
			weight=0.25;
		}

		if(JetTrueD0) {
			for(unsigned int id=0; id<TRUEDID->size(); ++id) {
				if(TMath::Abs(TRUEDID->at(id))!=421.) continue;
				if(TRUEDFROMB->at(id)) {
					trueD05->Fill(JetPT,weight);
				} else {
					trueD04->Fill(JetPT,weight);
				}
				break;
			}
		}
		if(JetTrueDSV) {
			if(SVN->size()>0) trueSV4->Fill(JetPT,weight);
		}
		if(JetTrueBSV) {
			if(SVN->size()>0) trueSV5->Fill(JetPT,weight);
		}
	}

	return true;
}

TH1D* unfold(TH1D* input, TString file, jetType type) {
	if(!input) return 0;

	TString nameStr="";

	switch(type) {
		case jetRecoD04:
			std::cout << "INFO : unfolding c2D0 results for file " << file << std::endl;
			nameStr="_c2D0";
			break;
		case jetRecoSV4:
			std::cout << "INFO : unfolding c2SV results for file " << file << std::endl;
			nameStr="_c2SV";
			break;
		case jetRecoD05:
			std::cout << "INFO : unfolding b2D0 results for file " << file << std::endl;
			nameStr="_b2D0";
			break;
		case jetRecoSV5:
			std::cout << "INFO : unfolding b2SV results for file " << file << std::endl;
			nameStr="_b2SV";
			break;
		default:
			return 0;
	}
	if(d0minpt!=5000. || d0maxpt!=-1) {
		nameStr+="_D0";
		nameStr+=d0minpt;
		nameStr+="-";
		nameStr+=d0maxpt;
	}

	//train unfolding on MC
	//use binning from data histogram
	RooUnfoldResponse* response = trainUnfold(input,type,file);

	if(!response) return 0;

	//plot response matrix
	TH1D* mcMeas  = dynamic_cast<TH1D*>(response->Hmeasured());
	TH1D* mcTrue  = dynamic_cast<TH1D*>(response->Htruth());
	TH2D* respMat = dynamic_cast<TH2D*>(response->Hresponse());

	mcTrue->SetLineStyle(kDashed);
	mcMeas->SetMaximum(1.1*TMath::Max(mcMeas->GetMaximum(),mcTrue->GetMaximum()));

	TCanvas c;
	respMat->Draw("colz");
	c.SaveAs(savedir+"/unfoldingResponse"+nameStr+".pdf");

	mcMeas->Draw();
	mcTrue->Draw("same");
	c.SaveAs(savedir+"/unfoldingMC"+nameStr+".pdf");

	//apply response
	RooUnfoldIds     unfold(response, input, 2);
	TH1D* ret = dynamic_cast<TH1D*>(unfold.Hreco());

	return ret;
}

int main(int argc, char** argv) {
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	//0-99 are MC closure tests, 100 is dijet tagged jets, 101 is jets with backwards SVs, 102 is J/psi tagged jets
	TString file="100";
	if(argc>1) file = argv[1];
	if(argc>2) savedir = argv[2];
	if(argc>3) d0minpt = atoi(argv[3]);//5000
	if(argc>4) d0maxpt = atoi(argv[4]);//-1

	gSystem->Exec("mkdir -p "+savedir);

	if(file.IsDec() && atoi(file)<100) {
		std::cout << "Running MC closure test " << file << " (" << atoi(file) << ")" << std::endl;
		makeTrainTestSamples(atoi(file));
	} else {
		if(file.BeginsWith("14") || file.BeginsWith("15") || file.BeginsWith("16") || file.BeginsWith("2") || file.BeginsWith("45")) {
			dataIsMC=true;
		}
		dataFile = "/eos/user/d/dcraik/jets-tuples-new-190210/for_yandex_data_new_"+file+".root";
	}

	unsigned int npt=3;
	double* binsPt  = new double[npt +1]{15000.,20000.,30000.,50000.};
	//double* binsPt  = new double[npt +1]{15000.,20000.,30000.,100000.};
	//unsigned int npt=4;
	//double* binsPt  = new double[npt +1]{10000.,15000.,20000.,30000.,100000.};

	//make input files
	if(useSimpleEff) {
		if(!addEffsSimple(file)) return 1;
	} else {
		if(!addEffs(file)) return 1;
	}
	if(!weightMC(file,jetRecoD04)) return 1;
	if(!weightMC(file,jetRecoD05)) return 1;
	if(!weightMC(file,jetRecoSV4)) return 1;
	if(!weightMC(file,jetRecoSV5)) return 1;

	//truth histograms for MC studies
	TH1D trueD04("trueD04","",npt,binsPt);
	TH1D trueD05("trueD05","",npt,binsPt);
	TH1D trueSV4("trueSV4","",npt,binsPt);
	TH1D trueSV5("trueSV5","",npt,binsPt);
	TH1D trueUnfldD04("trueUnfldD04","",npt,binsPt);
	TH1D trueUnfldD05("trueUnfldD05","",npt,binsPt);
	TH1D trueUnfldSV4("trueUnfldSV4","",npt,binsPt);
	TH1D trueUnfldSV5("trueUnfldSV5","",npt,binsPt);
	if(dataIsMC) {
		getTruth(&trueD04,&trueD05,&trueSV4,&trueSV5);
		getTruth(&trueUnfldD04,&trueUnfldD05,&trueUnfldSV4,&trueUnfldSV5,true);
		//for the real MC samples, test D0 efficiencies and return
		//for permutations, compare truth results to extracted yields
		if(!dataIsResampledMC) {
			if(useSimpleEff) {
				testEffsSimple(file,&trueD04);
			} else {
				testEffs(file,&trueD04);
			}
			return 0;//TODO
		}
	}

	//first do D0 fits for denominators
	TH1D recoD04("recoD04","",npt,binsPt);
	TH1D recoD05("recoD05","",npt,binsPt);
	if(!getD0Yields(&recoD04,&recoD05,file)) return 1;

	//do SV fits for numerator
	TH1D recoSV4("recoSV4","",npt,binsPt);
	TH1D recoSV5("recoSV5","",npt,binsPt);
	if(!getSVYields(&recoSV4,&recoSV5,file)) return 1;

	//do unfolding
	TH1D* unfoldedD04 = unfold(&recoD04, file, jetRecoD04);
	TH1D* unfoldedD05 = unfold(&recoD05, file, jetRecoD05);
	TH1D* unfoldedSV4 = unfold(&recoSV4, file, jetRecoSV4);
	TH1D* unfoldedSV5 = unfold(&recoSV5, file, jetRecoSV5);

	if(!unfoldedD04 || !unfoldedD05) return 1;
	if(!unfoldedSV4 || !unfoldedSV5) return 1;

	//print results
	double bfD0 = 0.0389;
	double errBFD0 = 0.0004;
	double ffc2D0 = 0.542;
	double errFFc2D0 = TMath::Sqrt(.024*.024 + .007*.007);
	double ffb2D0 = 0.587;
	double errFFb2D0 = TMath::Sqrt(.021*.021 + .008*.008);

	double bfffErr4 = TMath::Sqrt( (errFFc2D0/ffc2D0)*(errFFc2D0/ffc2D0) + (errBFD0/bfD0)*(errBFD0/bfD0) );
	double bfffErr5 = TMath::Sqrt( (errFFb2D0/ffb2D0)*(errFFb2D0/ffb2D0) + (errBFD0/bfD0)*(errBFD0/bfD0) );

	//D0 results
	std::cout << "D0 reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinError(i)/recoD04.GetBinContent(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinError(i)/recoD05.GetBinContent(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinContent(i) << " +/- " << recoD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinContent(i) << " +/- " << recoD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04->GetBinContent(i) << " +/- " << unfoldedD04->GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05->GetBinContent(i) << " +/- " << unfoldedD05->GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD04.GetBinContent(i) << " +/- " << trueD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD05.GetBinContent(i) << " +/- " << trueD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD04.GetBinContent(i) << " +/- " << trueUnfldD04.GetBinError(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD05.GetBinContent(i) << " +/- " << trueUnfldD05.GetBinError(i) << std::endl;
	std::cout << std::endl;

	//D0 results scaled
	std::cout << "jet reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << recoD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*recoD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << recoD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*recoD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << unfoldedD04->GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << unfoldedD05->GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << trueD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*trueD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << trueD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*trueD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << trueUnfldD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << trueUnfldD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	//SV results
	std::cout << "SV reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV4.GetBinContent(i) << " +/- " << recoSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV5.GetBinContent(i) << " +/- " << recoSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV4->GetBinContent(i) << " +/- " << unfoldedSV4->GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV5->GetBinContent(i) << " +/- " << unfoldedSV5->GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV4.GetBinContent(i) << " +/- " << trueSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV5.GetBinContent(i) << " +/- " << trueSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV4.GetBinContent(i) << " +/- " << trueUnfldSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV5.GetBinContent(i) << " +/- " << trueUnfldSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;

	std::ofstream fout("toysSVResults.log",std::ofstream::app);
	fout << file << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoSV4.GetBinContent(i) << "\t" << recoSV4.GetBinError(i) << "\t";
	for (unsigned int i=1; i<=npt; ++i) fout << recoSV5.GetBinContent(i) << "\t" << recoSV5.GetBinError(i) << "\t";
	fout << std::endl;
	fout.close();

	//ratios
	std::vector<double> ratioRec4;
	std::vector<double> ratioRec5;
	std::vector<double> ratioUnfld4;
	std::vector<double> ratioUnfld5;
	std::vector<double> ratioTrue4;
	std::vector<double> ratioTrue5;
	std::vector<double> ratioTrueUnfld4;
	std::vector<double> ratioTrueUnfld5;
	std::vector<double> errorRec4;
	std::vector<double> errorRec5;
	std::vector<double> errorUnfld4;
	std::vector<double> errorUnfld5;
	for (unsigned int i=1; i<=npt; ++i) {
		ratioRec4.push_back(recoSV4.GetBinContent(i) / (recoD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioRec5.push_back(recoSV5.GetBinContent(i) / (recoD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioUnfld4.push_back(unfoldedSV4->GetBinContent(i) / (unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioUnfld5.push_back(unfoldedSV5->GetBinContent(i) / (unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioTrue4.push_back(trueSV4.GetBinContent(i) / (trueD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioTrue5.push_back(trueSV5.GetBinContent(i) / (trueD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		ratioTrueUnfld4.push_back(trueUnfldSV4.GetBinContent(i) / (trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0)));
		ratioTrueUnfld5.push_back(trueUnfldSV5.GetBinContent(i) / (trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0)));
		errorRec4.push_back(ratioRec4[i-1]*TMath::Sqrt(TMath::Power(recoSV4.GetBinError(i)/recoSV4.GetBinContent(i),2)+TMath::Power(recoD04.GetBinError(i)/recoD04.GetBinContent(i),2)));
		errorRec5.push_back(ratioRec5[i-1]*TMath::Sqrt(TMath::Power(recoSV5.GetBinError(i)/recoSV5.GetBinContent(i),2)+TMath::Power(recoD05.GetBinError(i)/recoD05.GetBinContent(i),2)));
		errorUnfld4.push_back(ratioUnfld4[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV4->GetBinError(i)/unfoldedSV4->GetBinContent(i),2)+TMath::Power(unfoldedD04->GetBinError(i)/unfoldedD04->GetBinContent(i),2)));
		errorUnfld5.push_back(ratioUnfld5[i-1]*TMath::Sqrt(TMath::Power(unfoldedSV5->GetBinError(i)/unfoldedSV5->GetBinContent(i),2)+TMath::Power(unfoldedD05->GetBinError(i)/unfoldedD05->GetBinContent(i),2)));
	}
	std::cout << "ratios reco4/5/unfld4/5/true4/5/trueunfld4/5" << std::endl;
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioRec4[i] << "+/-" << errorRec4[i] << "+/-" << ratioRec4[i]*bfffErr4 << "\t";
	std::cout << std::endl;                                                              
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioRec5[i] << "+/-" << errorRec5[i] << "+/-" << ratioRec5[i]*bfffErr5 << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioUnfld4[i] << "+/-" << errorUnfld4[i] << "+/-" << ratioUnfld4[i]*bfffErr4 << "\t";
	std::cout << std::endl;
	for (unsigned int i=0; i<npt; ++i) std::cout << ratioUnfld5[i] << "+/-" << errorUnfld5[i] << "+/-" << ratioUnfld5[i]*bfffErr5 << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV4.GetBinContent(i) / (trueD04.GetBinContent(i) / (bfD0 * ffc2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueSV5.GetBinContent(i) / (trueD05.GetBinContent(i) / (bfD0 * ffb2D0)) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV4.GetBinContent(i) / (trueUnfldD04.GetBinContent(i) / (bfD0 * ffc2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << trueUnfldSV5.GetBinContent(i) / (trueUnfldD05.GetBinContent(i) / (bfD0 * ffb2D0)) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;

	//plots
	recoD04.SetMinimum(0.);
	recoD04.SetMaximum(1.1*TMath::Max(TMath::Max(recoD04.GetMaximum(),recoD05.GetMaximum()),TMath::Max(unfoldedD04->GetMaximum(),unfoldedD05->GetMaximum())));
	recoD04.SetLineColor(kBlue);
	recoD05.SetLineColor(kRed);
	unfoldedD04->SetLineColor(kBlue);
	unfoldedD05->SetLineColor(kRed);
	unfoldedD04->SetLineStyle(kDashed);
	unfoldedD05->SetLineStyle(kDashed);
	trueD04.SetLineColor(kBlue);
	trueD05.SetLineColor(kRed);
	trueD04.SetLineStyle(kDotted);
	trueD05.SetLineStyle(kDotted);

	recoSV4.SetMinimum(0.);
	recoSV4.SetMaximum(1.1*TMath::Max(TMath::Max(recoD04.GetMaximum(),recoD05.GetMaximum()),TMath::Max(unfoldedD04->GetMaximum(),unfoldedD05->GetMaximum())));
	recoSV4.SetLineColor(kBlue);
	recoSV5.SetLineColor(kRed);
	unfoldedSV4->SetLineColor(kBlue);
	unfoldedSV5->SetLineColor(kRed);
	unfoldedSV4->SetLineStyle(kDashed);
	unfoldedSV5->SetLineStyle(kDashed);
	trueSV4.SetLineColor(kBlue);
	trueSV5.SetLineColor(kRed);
	trueSV4.SetLineStyle(kDotted);
	trueSV5.SetLineStyle(kDotted);

	TCanvas c;
	recoD04.Draw();
	unfoldedD04->Draw("same");
	recoD05.Draw("same");
	unfoldedD05->Draw("same");
	if(dataIsMC) trueD04.Draw("same");
	if(dataIsMC) trueD05.Draw("same");
	c.SaveAs(savedir+"/d0Unfolding"+file+".pdf");

	recoSV4.Draw();
	unfoldedSV4->Draw("same");
	recoSV5.Draw("same");
	unfoldedSV5->Draw("same");
	if(dataIsMC) trueSV4.Draw("same");
	if(dataIsMC) trueSV5.Draw("same");
	c.SaveAs(savedir+"/svUnfolding"+file+".pdf");

	return 0;
}
