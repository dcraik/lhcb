#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
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

fitType whichFit(fitSBS1D);
double d0minpt(5000);
double d0maxpt(-1);

enum jetType {
	jetRecoD04,
	jetRecoD05,
	jetRecoSV4,
	jetRecoSV5
};

// function to make plot of 1D projection of the fit
void plot(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title) {
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

	c1.SaveAs("plots/"+name+"_fit.pdf");
	c1.SaveAs("plots/"+name+"_fit.png");

	delete plot;
}

// function to print parameters to a file
void print(TString file, const RooArgList& params) {
	FILE * pFile = fopen(file.Data(), "w");

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

RooUnfoldResponse* trainUnfold(TH1* binning, jetType type, TString file) {
	TFile* f(0);
	TString nameStr;
	switch(type) {
		case jetRecoD04:
			nameStr="_c2D0";
			break;
		case jetRecoSV4:
			nameStr="_c2SV";
			break;
		case jetRecoD05:
			nameStr="_b2D0";
			break;
		case jetRecoSV5:
			nameStr="_b2SV";
			break;
		default:
			return 0;
	}
	f = TFile::Open("D0MCjets"+file+nameStr+".root");
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


bool addEffs(TString file) {
	TFile* f0 = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-181001/for_yandex_data_new_"+file+".root");

	TFile* f2 = TFile::Open("../efficiencies/D0Effs_2XX_16x8bins_splitEffs.root");
	TFile* f3 = TFile::Open("../efficiencies/D0AccEffNew.root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	if(!t0) return false;

	TH2D* hacc = dynamic_cast<TH2D*>(f3->Get("eff"));
	TH2D* hrec = dynamic_cast<TH2D*>(f3->Get("reff"));

	TH2D* hpid  = dynamic_cast<TH2D*>(f2->Get("pideffD045"));
	TH2D* hsel4 = dynamic_cast<TH2D*>(f2->Get("seleffD04"));
	TH2D* hsel5 = dynamic_cast<TH2D*>(f2->Get("seleffD05"));
	//TH2D* haccB = dynamic_cast<TH2D*>(f2->Get("acceffD045"));
	//TH2D* hrecB = dynamic_cast<TH2D*>(f2->Get("receffD045"));

	if(!hacc || !hrec || !hpid || !hsel4 || !hsel5) return false;

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

	unsigned int nentries0 = t0->GetEntries();

	TFile* fout = TFile::Open("D0jets"+file+".root","RECREATE");
	TTree* tout = new TTree("T","");

	double D0M(0.), D0Pt(0.), D0Eta(0.), D0LogIPChi2(0.), weight4(0.), weight5(0.);

	tout->Branch("JetPT",        &JetPT);
	tout->Branch("JetEta",       &JetEta);
	tout->Branch("D0M",          &D0M);
	tout->Branch("D0PT",         &D0Pt);
	tout->Branch("D0Eta",        &D0Eta);
	tout->Branch("D0LogIPChi2",  &D0LogIPChi2);
	tout->Branch("weight4",      &weight4);
	tout->Branch("weight5",      &weight5);

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
			if(vD0KPNNK->at(s)<0.2 || vD0PIPNNPI->at(s)<0.1) continue;

			D0M         = vD0M->at(s);
			D0LogIPChi2 = TMath::Log(vD0IPCHI2->at(s));
			D0Pt        = D0P.Pt();
			D0Eta       = D0P.Eta();

			//get efficiency corrections
			double effacc(0.), effrec(0.), effsel4(0.), effsel5(0.), effpid(0.);

			if(D0Pt>=100000.) D0Pt=99999.;
			effacc =      hacc ->GetBinContent(hacc ->FindBin(D0Pt/1000.,D0Eta));
			effrec =      hrec ->GetBinContent(hrec ->FindBin(D0Pt/1000.,D0Eta));
			effsel4=      hsel4->GetBinContent(hsel4->FindBin(D0Pt      ,D0Eta));
			effsel5=      hsel5->GetBinContent(hsel5->FindBin(D0Pt      ,D0Eta));
			effpid =      hpid ->GetBinContent(hpid ->FindBin(D0Pt      ,D0Eta));

			double eff4=effacc*effrec*effsel4*effpid;
			double eff5=effacc*effrec*effsel5*effpid;

			if(eff4<0.01 && eff5<0.01) {
				std::cout << D0Pt << "\t" << D0Eta << "\t" << effacc << "\t" << effrec << "\t" << effsel4 << "\t" << effsel5 << "\t" << effpid << std::endl;
				continue;
			}

			weight4 = 1./eff4;
			weight5 = 1./eff5;

			tout->Fill();
			break;//only keep one D0 candidate per entry
		}
	}

	tout->AutoSave();
	fout->Close();

	return true;
}

bool weightMC(TString file, jetType type) {
	TString nameStr="";

	TFile* fin(0);
	switch(type) {
		case jetRecoD04:
			nameStr="_c2D0";
			fin = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_14X.root");
			break;
		case jetRecoSV4:
			nameStr="_c2SV";
			fin = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_14X.root");
			break;
		case jetRecoD05:
			nameStr="_b2D0";
			fin = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_15X.root");
			break;
		case jetRecoSV5:
			nameStr="_b2SV";
			fin = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_15X.root");
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
		c.SaveAs("D0fPtReweighting"+file+nameStr+".pdf");

		fPtDWeights->Draw();
		c.SaveAs("D0fPtWeights"+file+nameStr+".pdf");
	}
	//end f(pT) weights setup

	//now make the weighted output file
	TFile* fout = TFile::Open("D0MCjets"+file+nameStr+".root","RECREATE");
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
		if(type==jetRecoD04||type==jetRecoD05)
			weight *= fPtDWeights->GetBinContent(fPtDWeights->FindBin(D0PT/JetPT));

		tout->Fill();
	}

	tout->SetDirectory(fout);
	tout->Write();
	tout->AutoSave();
	fout->Close();

	return true;
}

bool fitD0_1D1D(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
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
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));

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
	//	promptMean.setConstant();
	//	promptWidth.setConstant();
		promptAsym.setConstant();
		promptRhoL.setConstant();
		promptRhoR.setConstant();
	//	displacedMean.setConstant();
	//	displacedWidth.setConstant();
	}

	RooFitResult * result2 = data_pdf2_full.fitTo( *dsPeak, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE), RooFit::Constrain(RooArgSet(bkgYield)));

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

	plot(DM, massVal-massScaleFact*80., massVal+massScaleFact*80., 40, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_M_"+ptStr, "m");
	plot(DLOGIPCHI2, -5., 15., 20, dsPeak, data_pdf2, sig_pdfs2, bkg_pdfs2, dType+"_"+file+"_IPChi2_"+ptStr, "log(IP#chi^{2})");

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
	print(dType+"_"+file+"_params_"+ptStr+".dat",params);

	return true;
}

bool fitD0_SBS(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
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
	c.SaveAs(dType+"_"+file+"_IPChi2_SBS_"+ptStr+".pdf");

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

	hsub.Add(&hpeak,&hside,1.,-1.);

	TCanvas c;
	hpeak.Draw();
	hsub.SetLineColor(kGreen+2);
	hsub.Draw("same");
	hside.SetLineColor(kRed);
	hside.Draw("same");
	c.SaveAs(dType+"_"+file+"_IPChi2_SBS_"+ptStr+".pdf");

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

	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));

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

	plot(DLOGIPCHI2, -5., 15., 20, &dh, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_IPChi2_SBS1D_"+ptStr, "log(IP#chi^{2})");

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
	print(dType+"_"+file+"_params_SBS1D_"+ptStr+".dat",params);

	return true;
}

bool fitD0_2D(int flavour, double minpt, double maxpt, double& yield, double& error, TString file) {
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
	/*RooFitResult * result =*/ data_pdf.fitTo(*ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::SumW2Error(kTRUE));

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

	plot(DM, massVal-massScaleFact*80., massVal+massScaleFact*80., 40, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_M_2D_"+ptStr, "m");
	plot(DLOGIPCHI2, -5., 15., 20, ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"_"+file+"_IPChi2_2D_"+ptStr, "log(IP#chi^{2})");

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
	print(dType+"_"+file+"_params2D_"+ptStr+".dat",params);

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
		hist4->SetBinContent(i,yield);
		hist4->SetBinError(  i,error);
	}

	for(int i=1; i<=hist5->GetNbinsX(); ++i) {
		if(!fitD0(5,hist5->GetBinLowEdge(i),hist5->GetBinLowEdge(i+1),yield,error,file)) return false;
		hist5->SetBinContent(i,yield);
		hist5->SetBinError(  i,error);
	}

	return true;
}

//function to fit features for a single sample
TString fitSV(double& NC, double& eC, double& NB, double& eB, TString sample="100", double minPT=20000, double maxPT=30000){

	TString ptCutStr;
	ptCutStr = "JetPT>"; ptCutStr+=minPT; ptCutStr+=" && JetPT<"; ptCutStr+=maxPT;

	TFile* f0 = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_101.root");
	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));

	TFile* f4 = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_14X.root");
	TTree* t4 = dynamic_cast<TTree*>(f4->Get("T"));

	TFile* f5 = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180914/for_yandex_data_new_15X.root");
	TTree* t5 = dynamic_cast<TTree*>(f5->Get("T"));

	TFile* fd = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-181001/for_yandex_data_new_"+sample+".root");
	TTree* td = dynamic_cast<TTree*>(fd->Get("T"));


	if(!t0 || !t4 || !t5 || !td) return "FAIL";

	TH1D svmcorsvn_0("svmcorsvn_0","",28,0.5,28.5);
	TH1D svmcorsvn_4("svmcorsvn_4","",28,0.5,28.5);
	TH1D svmcorsvn_5("svmcorsvn_5","",28,0.5,28.5);
	TH1D svmcorsvn_d("svmcorsvn_d","",28,0.5,28.5);

	svmcorsvn_4.Sumw2();
	svmcorsvn_5.Sumw2();
	
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

	//get MC/backwards efficiencies for estimating total yields
	int eff0denom(0.), eff4denom(0.), eff5denom(0.);
	int eff0num(0.), eff4num(0.), eff5num(0.);
	double weight(1.);

	unsigned int npt=4;
	double* ptInputWeightBins  = new double[npt +1]{10000.,15000.,20000.,50000.,100000.};

	TH1D* mcWeights4 = new TH1D("mcWeights4","",npt,ptInputWeightBins);
	TH1D* mcWeights5 = new TH1D("mcWeights5","",npt,ptInputWeightBins);
	mcWeights4->SetBinContent(1,1.);
	mcWeights4->SetBinContent(2,0.3);
	mcWeights4->SetBinContent(3,0.12);
	mcWeights4->SetBinContent(4,0.008);
	mcWeights5->SetBinContent(1,1.);
	mcWeights5->SetBinContent(2,0.2);
	mcWeights5->SetBinContent(3,0.08);
	mcWeights5->SetBinContent(4,0.006);

	for(int i=0; i<t0->GetEntries(); ++i) {
		t0->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff0denom;
		if(SVN->size()<1) continue;
		++eff0num;
		if(SVN->at(0)>1.5 && SVN->at(0)<10.5) svmcorsvn_0.Fill(SVN->at(0)+18.);
		else if(SVN->at(0)>10.5) svmcorsvn_0.Fill(28.);
		if(SVMCor->at(0)>500. && SVMCor->at(0)<10000.) svmcorsvn_0.Fill(SVMCor->at(0)/500. -0.5);
		else if(SVMCor->at(0)>10000.) svmcorsvn_0.Fill(19.);
	}
	std::cout << svmcorsvn_0.Integral(1,19) << "\t" << svmcorsvn_0.Integral(20,28) << std::endl;

	for(int i=0; i<t4->GetEntries(); ++i) {
		t4->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff4denom;
		if(SVN->size()<1) continue;
		++eff4num;
		if(JetTruePT>=100000.) JetTruePT=99999.;
		weight = mcWeights4->GetBinContent(mcWeights4->FindBin(JetTruePT));
		if(SVN->at(0)>1.5 && SVN->at(0)<10.5) svmcorsvn_4.Fill(SVN->at(0)+18.,weight);
		else if(SVN->at(0)>10.5) svmcorsvn_4.Fill(28.,weight);
		if(SVMCor->at(0)>500. && SVMCor->at(0)<10000.) svmcorsvn_4.Fill(SVMCor->at(0)/500. -0.5,weight);
		else if(SVMCor->at(0)>10000.) svmcorsvn_4.Fill(19.,weight);
	}
	std::cout << svmcorsvn_4.Integral(1,19) << "\t" << svmcorsvn_4.Integral(20,28) << std::endl;

	for(int i=0; i<t5->GetEntries(); ++i) {
		t5->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		++eff5denom;
		if(SVN->size()<1) continue;
		++eff5num;
		if(JetTruePT>=100000.) JetTruePT=99999.;
		weight = mcWeights5->GetBinContent(mcWeights5->FindBin(JetTruePT));
		if(SVN->at(0)>1.5 && SVN->at(0)<10.5) svmcorsvn_5.Fill(SVN->at(0)+18.,weight);
		else if(SVN->at(0)>10.5) svmcorsvn_5.Fill(28.,weight);
		if(SVMCor->at(0)>0. && SVMCor->at(0)<10000.) svmcorsvn_5.Fill(SVMCor->at(0)/500. -0.5,weight);
		else if(SVMCor->at(0)>10000.) svmcorsvn_5.Fill(19.,weight);
	}
	std::cout << svmcorsvn_5.Integral(1,19) << "\t" << svmcorsvn_5.Integral(20,28) << std::endl;

	for(int i=0; i<td->GetEntries(); ++i) {
		td->GetEntry(i);
		if(JetPT<minPT || JetPT>maxPT) continue;
		if(SVN->size()<1) continue;
		if(SVN->at(0)>1.5 && SVN->at(0)<10.5) svmcorsvn_d.Fill(SVN->at(0)+18.);
		else if(SVN->at(0)>10.5) svmcorsvn_d.Fill(28.);
		if(SVMCor->at(0)>0. && SVMCor->at(0)<10000.) svmcorsvn_d.Fill(SVMCor->at(0)/500. -0.5);
		else if(SVMCor->at(0)>10000.) svmcorsvn_d.Fill(19.);
	}
	std::cout << svmcorsvn_d.Integral(1,19) << "\t" << svmcorsvn_d.Integral(20,28) << std::endl;

	// -- variables from datasets
	RooRealVar SVComb(  "SVComb",  "SVComb",  0.5,  28.5,  ""); 
	SVComb.setBins(28);

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
	/*RooFitResult * result =*/ data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfQ" );
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	std::vector<std::string> bkg_pdfs;

	TString plotName = "SVComb_"+sample+"_"; plotName+=minPT; plotName+="-"; plotName+=maxPT;
	plot(SVComb, 0.5, 28.5, 28, &dh, data_pdf, sig_pdfs, bkg_pdfs, plotName, "M_{cor}, N_{trk}");

	//print parameters
	RooArgList params;
	params.add(yieldB);
	params.add(yieldC);
	params.add(yieldQ);
	TString paramsName = "params_comb_"+sample+"_"; paramsName+=minPT; paramsName+="-"; paramsName+=maxPT; paramsName+=".dat";
	print(paramsName,params);

	double NQ, Ntot;
	double eQ;

	Ntot = yieldB.getValV()+yieldC.getValV()+yieldQ.getValV();
	NB = yieldB.getValV()/2.;
	NC = yieldC.getValV()/2.;
	NQ = yieldQ.getValV()/2.;
	eB = yieldB.getError()/2.;
	eC = yieldC.getError()/2.;
	eQ = yieldQ.getError()/2.;
	Ntot/=2.;
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
		hist4->SetBinContent(i,nC);
		hist4->SetBinError(  i,eC);
		hist5->SetBinContent(i,nB);
		hist5->SetBinError(  i,eB);
	}

	return true;
}

TH1D* unfold(TH1D* input, TString file, jetType type) {
	if(!input) return 0;

	TString nameStr="";

	switch(type) {
		case jetRecoD04:
			nameStr="_c2D0";
			break;
		case jetRecoSV4:
			nameStr="_c2SV";
			break;
		case jetRecoD05:
			nameStr="_b2D0";
			break;
		case jetRecoSV5:
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
	c.SaveAs("unfoldingResponse"+nameStr+".pdf");

	mcMeas->Draw();
	mcTrue->Draw("same");
	c.SaveAs("unfoldingMC"+nameStr+".pdf");

	//apply response
	RooUnfoldIds     unfold(response, input, 2);
	TH1D* ret = dynamic_cast<TH1D*>(unfold.Hreco());

	return ret;
}

int main(int argc, char** argv) {
	gStyle->SetOptStat(0);

	TString file="100";
	if(argc>1) file = argv[1];
	if(argc>2) d0minpt = atoi(argv[2]);//5000
	if(argc>3) d0maxpt = atoi(argv[3]);//-1

	unsigned int npt=3;
	double* binsPt  = new double[npt +1]{15000.,20000.,30000.,100000.};
	//unsigned int npt=4;
	//double* binsPt  = new double[npt +1]{10000.,15000.,20000.,30000.,100000.};
	
	//make input files
	if(!addEffs(file)) return 1;
	if(!weightMC(file,jetRecoD04)) return 1;
	if(!weightMC(file,jetRecoD05)) return 1;
	if(!weightMC(file,jetRecoSV4)) return 1;
	if(!weightMC(file,jetRecoSV5)) return 1;

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
	std::cout << "D0" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinError(i)/recoD04.GetBinContent(i) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinError(i)/recoD05.GetBinContent(i) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD04.GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << recoD04.GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*recoD04.GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoD05.GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << recoD05.GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*recoD05.GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;
	std::cout << std::endl;

	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0) << " +/- " << unfoldedD04->GetBinError(i) / (bfD0 * ffc2D0) << " +/- " << bfffErr4*unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0) << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0) << " +/- " << unfoldedD05->GetBinError(i) / (bfD0 * ffb2D0) << " +/- " << bfffErr5*unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0) << std::endl;

	//SV results
	std::cout << "SV" << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV4.GetBinContent(i) << " +/- " << recoSV4.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV4->GetBinContent(i) << " +/- " << unfoldedSV4->GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV5.GetBinContent(i) << " +/- " << recoSV5.GetBinError(i) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV5->GetBinContent(i) << " +/- " << unfoldedSV5->GetBinError(i) << "\t";
	std::cout << std::endl;
	std::cout << std::endl;

	//ratios
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV4.GetBinContent(i) / (recoD04.GetBinContent(i) / (bfD0 * ffc2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV4->GetBinContent(i) / (unfoldedD04->GetBinContent(i) / (bfD0 * ffc2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << recoSV5.GetBinContent(i) / (recoD05.GetBinContent(i) / (bfD0 * ffb2D0)) << "\t";
	std::cout << std::endl;
	for (unsigned int i=1; i<=npt; ++i) std::cout << unfoldedSV5->GetBinContent(i) / (unfoldedD05->GetBinContent(i) / (bfD0 * ffb2D0)) << "\t";
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

	TCanvas c;
	recoD04.Draw();
	unfoldedD04->Draw("same");
	recoD05.Draw("same");
	unfoldedD05->Draw("same");
	c.SaveAs("dataUnfolding"+file+".pdf");

	return 0;
}
