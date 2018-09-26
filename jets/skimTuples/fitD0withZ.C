#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTree.h"

#include "RooAbsDataStore.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include "RooStats/SPlot.h"

#include <vector>
#include <fstream>

#include "RooPromptShape.h"

// function to make plot of 1D projection of the fit
void plot(RooRealVar& var, double min, double max, RooDataSet& dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title) {
	TCanvas c1;
	c1.SetBottomMargin(0.19);
	RooPlot* plot = var.frame(min, max);
	dh.plotOn( plot );
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
	dh.plotOn( plot );
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

//function to fit features for a single sample
void fit(){

	// -- variables from datasets
	RooRealVar Z0M(  "Z0M",  "Z0M",  40000.,  120000.,  "MeV/#it{c}^{2}"); 
	Z0M.setBins(40);
	RooRealVar D0M(  "D0M",  "D0M",  1784.,  1944.,  "MeV/#it{c}^{2}"); 
	D0M.setBins(30);
	RooRealVar D0LOGIPCHI2(  "D0LOGIPCHI2",  "D0LOGIPCHI2",  -5.,  15.,  ""); 
	D0LOGIPCHI2.setBins(20);

	// -- add all feature observables to dataset
	RooRealVar DE(  "D0E",  "D0E",  0.,   2e6); 
	RooRealVar DPX( "D0PX", "D0PX", -1e5, 1e5); 
	RooRealVar DPY( "D0PY", "D0PY", -1e5, 1e5); 
	RooRealVar DPZ( "D0PZ", "D0PZ", 0.,   2e6); 

	RooArgSet obs;
	obs.add(Z0M);
	obs.add(D0M);
	obs.add(D0LOGIPCHI2);
	obs.add(DE);
	obs.add(DPX);
	obs.add(DPY);
	obs.add(DPZ);

	// -- load dataset for data
	TFile* df = TFile::Open("/tmp/dcraik/ZjetSkimmed.root");
	TTree* dt = dynamic_cast<TTree*>(df->Get("T"));
	RooDataSet ds("ds","ds", obs, RooFit::Import(*dt));
	RooDataSet* dsPeak = dynamic_cast<RooDataSet*>(ds.reduce("D0M>1844 && D0M<1884"));
	RooDataSet* dsSB   = dynamic_cast<RooDataSet*>(ds.reduce("D0M<1804 || D0M>1924"));

	//Z0M model
	RooRealVar zMass("zMass","mass mean",91000.,80000.,100000.);
	RooRealVar zWidth("zWidth","mass width",2300.,1000.,6000.);
	RooGaussian sigZMass("sigZMass","",Z0M,zMass,zWidth);

	RooRealVar p0Z("p0Z","p0Z",0., -0.0001, 0.0001);
	RooPolynomial bkgZMass("bkgZMass","",Z0M, RooArgList(p0Z));

	// -- yields
	RooRealVar zSigYield(  "zSigYield",  "yield Z",    50,    0.0,     100);
	RooRealVar zBkgYield(  "zBkgYield",  "yield comb Z", 20,    0.0,     100);

	// -- total PDF
	RooAddPdf Zdata_pdf( "Zdata_pdf",  "Zdata_pdf", RooArgList(sigZMass,bkgZMass), RooArgList(zSigYield,zBkgYield) );

	// -- fit model pdf to the dataset ----------------------------------------------
	RooFitResult * Zresult = Zdata_pdf.fitTo( ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

        Z0M.setRange("zsignal",91000.-2*2300.,91000.+2*2300.);
        Z0M.setRange("zsideLo",91000.-6*2300.,91000.-4*2300.);
        Z0M.setRange("zsideHi",91000.+4*2300.,91000.+6*2300.);
        Z0M.setRange("zfull",40000.,120000.);

        double fzsig_1 = sigZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zsignal"))->getVal();
        double fzsig_2 = sigZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zsideLo"))->getVal();
        double fzsig_3 = sigZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zsideHi"))->getVal();
        double fzsig_0 = sigZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zfull"))->getVal();

        double fzbkg_1 = bkgZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zsignal"))->getVal();
        double fzbkg_2 = bkgZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zsideLo"))->getVal();
        double fzbkg_3 = bkgZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zsideHi"))->getVal();
        double fzbkg_0 = bkgZMass.createIntegral(RooArgSet(Z0M),RooFit::NormSet(Z0M),RooFit::Range("zfull"))->getVal();

        std::cout << "Z yields in mass windows" << std::endl;
        std::cout << zSigYield.getVal()*fzsig_1/fzsig_0 << "\t" << zBkgYield.getVal()*fzbkg_1/fzbkg_0 << std::endl;
        std::cout << zSigYield.getVal()*(fzsig_2+fzsig_3)/fzsig_0 << "\t" << zBkgYield.getVal()*(fzbkg_2+fzbkg_3)/fzbkg_0 << std::endl;

	//D0M model
	RooRealVar dMass("dMass","mass mean",1864.);//,1844.,1884.);
	RooRealVar dWidth("dWidth","mass width",10.);//,5.,20.);
	RooGaussian sigMass("sigMass","",D0M,dMass,dWidth);

	RooRealVar p0("p0","p0",0., -0.001, 0.001);
	RooPolynomial bkgMass("bkgMass","",D0M, RooArgList(p0));

	// -- yields
	RooRealVar sigYield(  "sigYield",  "yield D",    20,    0.0,     100);
	RooRealVar bkgYield(  "bkgYield",  "yield comb", 40,    0.0,     100);

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(sigMass,bkgMass), RooArgList(sigYield,bkgYield) );

	// -- fit model pdf to the dataset ----------------------------------------------
	RooFitResult * result = data_pdf.fitTo( ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	sigYield.setConstant();
	bkgYield.setConstant();
	dMass.setConstant();
	dWidth.setConstant();
	p0.setConstant();

        D0M.setRange("signal",1844.,1884.);
        D0M.setRange("sideLo",1784.,1804.);
        D0M.setRange("sideHi",1924.,1944.);
        D0M.setRange("full",1784.,1944.);

        double fsig_1 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("signal"))->getVal();
        double fsig_2 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideLo"))->getVal();
        double fsig_3 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideHi"))->getVal();
        double fsig_0 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("full"))->getVal();

        double fbkg_1 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("signal"))->getVal();
        double fbkg_2 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideLo"))->getVal();
        double fbkg_3 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideHi"))->getVal();
        double fbkg_0 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("full"))->getVal();

        std::cout << "yields in mass windows" << std::endl;
        std::cout << sigYield.getVal()*fsig_1/fsig_0 << "\t" << bkgYield.getVal()*fbkg_1/fbkg_0 << std::endl;
        std::cout << sigYield.getVal()*(fsig_2+fsig_3)/fsig_0 << "\t" << bkgYield.getVal()*(fbkg_2+fbkg_3)/fbkg_0 << std::endl;
        //std::cout << sigYield.getVal()*fsig_2/fsig_0 << "\t" << bkgYield.getVal()*fbkg_2/fbkg_0 << std::endl;
        //std::cout << sigYield.getVal()*fsig_3/fsig_0 << "\t" << bkgYield.getVal()*fbkg_3/fbkg_0 << std::endl;

	RooRealVar promptMean("promptMean","mean prompt",0.,-1.,3.);
	//RooRealVar promptWidth("promptWidth","width prompt",1.,0.5,5.);
	//RooGaussian promptLOGIPCHI2("promptLOGIPCHI2","",D0LOGIPCHI2,promptMean,promptWidth);
	RooRealVar promptWidth("promptWidth","width prompt",1.,0.5,5.);
	RooRealVar promptAsym("promptAsym","asym. prompt",0.,-1.,1.);
	RooRealVar promptRhoL("promptRhoL","exp L prompt",1.3,0.5,2.);
	RooRealVar promptRhoR("promptRhoR","exp R prompt",1.3,0.5,2.);
	RooPromptShape promptLOGIPCHI2("promptLOGIPCHI2","",D0LOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	RooRealVar displacedMean("displacedMean","mean displ.",5.,3.,12.);
	RooRealVar displacedWidth("displacedWidth","width displ.",3.,1.,8.);
	RooGaussian displacedLOGIPCHI2("displacedLOGIPCHI2","",D0LOGIPCHI2,displacedMean,displacedWidth);

	RooKeysPdf bkgLOGIPCHI2("bkgLOGIPCHI2","",D0LOGIPCHI2,*dsSB);

	RooRealVar fSigInPeak("fSigInPeak","",fsig_1/fsig_0);
	RooRealVar fBkgInPeak("fBkgInPeak","",fbkg_1/fbkg_0);
	RooRealVar fPrompt("fPrompt","frac. prompt",0.5,0.,1.);
	RooFormulaVar promptYield("promptYield","","@0*@1*@2",RooArgList(fPrompt,fSigInPeak,sigYield));
	RooFormulaVar displacedYield("displacedYield","","(1.0-@0)*@1*@2",RooArgList(fPrompt,fSigInPeak,sigYield));
	RooFormulaVar bkgInPeakYield("bkgInPeakYield","","@0*@1",RooArgList(fBkgInPeak,bkgYield));

	RooAddPdf data_pdf2( "data_pdf2",  "data_pdf2", RooArgList(promptLOGIPCHI2,displacedLOGIPCHI2,bkgLOGIPCHI2), RooArgList(promptYield,displacedYield,bkgInPeakYield) );

	RooFitResult * result2 = data_pdf2.fitTo( *dsPeak, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	promptMean.setConstant();
	promptWidth.setConstant();
	promptAsym.setConstant();
	promptRhoL.setConstant();
	promptRhoR.setConstant();
	displacedMean.setConstant();
	displacedWidth.setConstant();

	//make plots
	std::vector<std::string> sig_pdfsZ;
	sig_pdfsZ.push_back( "sigZMass" );
	std::vector<std::string> bkg_pdfsZ;
	bkg_pdfsZ.push_back( "bkgZMass" );

	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "sigMass" );
	std::vector<std::string> bkg_pdfs;
	bkg_pdfs.push_back( "bkgMass" );

	std::vector<std::string> sig_pdfs2;
	sig_pdfs2.push_back( "promptLOGIPCHI2" );
	sig_pdfs2.push_back( "displacedLOGIPCHI2" );
	std::vector<std::string> bkg_pdfs2;
	bkg_pdfs2.push_back( "bkgLOGIPCHI2" );

	plot(Z0M,  40000., 120000., ds, Zdata_pdf, sig_pdfsZ, bkg_pdfsZ, "Z0M_Zjet", "m_{#mu#mu}");
	plot(D0M, 1784., 1944., ds, data_pdf, sig_pdfs, bkg_pdfs, "D0M_Zjet", "m_{K#pi}");
	plot(D0LOGIPCHI2, -5., 15., *dsPeak, data_pdf2, sig_pdfs2, bkg_pdfs2, "D0IPChi2_Zjet", "log(IP#chi^{2})");

	//print parameters
	RooArgList params;
	params.add(zSigYield);
	params.add(zBkgYield);
	params.add(zMass);
	params.add(zWidth);
	params.add(sigYield);
	params.add(bkgYield);
	params.add(fPrompt);
	params.add(dMass);
	params.add(dWidth);
	params.add(p0);
	params.add(promptMean);
	params.add(promptWidth);
	params.add(promptAsym);
	params.add(promptRhoL);
	params.add(promptRhoR);
	params.add(displacedMean);
	params.add(displacedWidth);
	print("ZD0params_Zjet.dat",params);

	//Start SPlot of file0
	//SPlot only works with RooRealVar objects so need to define a new PDF
	RooRealVar promptYieldSP(     "promptYieldSP",    "yield prompt", promptYield.getVal(),    0.0,     1e5);
	RooRealVar displacedYieldSP(  "displacedYieldSP", "yield displ.", displacedYield.getVal(), 0.0,     1e5);
	RooRealVar bkgYieldSP(        "bkgYieldSP",       "yield comb",   bkgInPeakYield.getVal(), 0.0,     1e5);
	RooAddPdf data_pdf2SP( "data_pdf2SP",  "data_pdf2SP", RooArgList(promptLOGIPCHI2,displacedLOGIPCHI2,bkgLOGIPCHI2), RooArgList(promptYieldSP,displacedYieldSP,bkgYieldSP) );
        RooStats::SPlot* sData0 = new RooStats::SPlot("sData","An SPlot",
                        *dsPeak, &data_pdf2SP, RooArgList(promptYieldSP, displacedYieldSP, bkgYieldSP) );
        RooDataSet * dsSW0 = new RooDataSet("dsSW0","dsSW0",dsPeak,*(dsPeak->get()),0);//,"promptYieldSP_sw");

        TTree * tree_data0 = (TTree*)dsSW0->store()->tree();
        TFile * newFile0 = TFile::Open("Z_D0_sWeights.root","RECREATE");
        tree_data0->SetName("T");
        tree_data0->Write();
        newFile0->Save();
        newFile0->Close();
	//Finish SPlot
}

int main() {
	fit();
}
