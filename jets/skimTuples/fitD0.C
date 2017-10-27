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
	RooRealVar D0M(  "D0M",  "D0M",  1784.,  1944.,  "MeV/#it{c}^{2}"); 
	D0M.setBins(30);
	RooRealVar D0LOGIPCHI2(  "D0LOGIPCHI2",  "D0LOGIPCHI2",  -5.,  15.,  ""); 
	D0LOGIPCHI2.setBins(20);

	RooRealVar dMass("dMass","mass mean",1864.,1844.,1884.);
	RooRealVar dWidth("dWidth","mass width",20.,5.,50.);
	RooGaussian sigMass("sigMass","",D0M,dMass,dWidth);

	RooRealVar p0("p0","p0",0., -0.1, 0.1);
	RooPolynomial bkgMass("bkgMass","",D0M, RooArgList(p0));

	// -- yields
	RooRealVar sigYield(  "sigYield",  "yield D",    750,    0.0,     1e5);
	RooRealVar bkgYield(  "bkgYield",  "yield comb", 750,    0.0,     1e5);

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(sigMass,bkgMass), RooArgList(sigYield,bkgYield) );

	// -- add all feature observables to dataset
	RooRealVar DE(  "D0E",  "D0E",  0.,   2e6); 
	RooRealVar DPX( "D0PX", "D0PX", -1e5, 1e5); 
	RooRealVar DPY( "D0PY", "D0PY", -1e5, 1e5); 
	RooRealVar DPZ( "D0PZ", "D0PZ", 0.,   2e6); 

	RooArgSet obs;
	obs.add(D0M);
	obs.add(D0LOGIPCHI2);
	obs.add(DE);
	obs.add(DPX);
	obs.add(DPY);
	obs.add(DPZ);

	// -- load dataset for data
	TFile* df = TFile::Open("for_yandex_data_SV_0tag_new170911_DsOnly_forFit.root");
	TTree* dt = dynamic_cast<TTree*>(df->Get("T"));
	RooDataSet ds("ds","ds", obs, RooFit::Import(*dt));
	RooDataSet* dsPeak = dynamic_cast<RooDataSet*>(ds.reduce("D0M>1844 && D0M<1884"));
	RooDataSet* dsSB   = dynamic_cast<RooDataSet*>(ds.reduce("D0M<1804 || D0M>1924"));

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
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "sigMass" );
	std::vector<std::string> bkg_pdfs;
	bkg_pdfs.push_back( "bkgMass" );

	std::vector<std::string> sig_pdfs2;
	sig_pdfs2.push_back( "promptLOGIPCHI2" );
	sig_pdfs2.push_back( "displacedLOGIPCHI2" );
	std::vector<std::string> bkg_pdfs2;
	bkg_pdfs2.push_back( "bkgLOGIPCHI2" );

	plot(D0M, 1784., 1944., ds, data_pdf, sig_pdfs, bkg_pdfs, "D0M_full", "m_{K#pi}");
	plot(D0LOGIPCHI2, -5., 15., *dsPeak, data_pdf2, sig_pdfs2, bkg_pdfs2, "D0IPChi2_full", "log(IP#chi^{2})");

	//print parameters
	RooArgList params;
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
	print("D0params_full.dat",params);

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
        TFile * newFile0 = TFile::Open("data0_D0_sWeights.root","RECREATE");
        tree_data0->SetName("T");
        tree_data0->Write();
        newFile0->Save();
        newFile0->Close();
	//Finish SPlot

	TFile* df4 = TFile::Open("for_yandex_data_SV_4tag_new170911_DsOnly_forFit.root");
	TTree* dt4 = dynamic_cast<TTree*>(df4->Get("T"));
	RooDataSet ds4("ds4","ds4", obs, RooFit::Import(*dt4));
	RooDataSet* ds4Peak = dynamic_cast<RooDataSet*>(ds4.reduce("D0M>1844 && D0M<1884"));

	sigYield.setConstant(false);
	bkgYield.setConstant(false);
	dMass.setConstant(false);
	dWidth.setConstant(false);
	p0.setConstant(false);

	RooFitResult * result3 = data_pdf.fitTo( ds4, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	sigYield.setConstant();
	bkgYield.setConstant();
	p0.setConstant();

        fsig_1 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("signal"))->getVal();
        fsig_2 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideLo"))->getVal();
        fsig_3 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideHi"))->getVal();
        fsig_0 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("full"))->getVal();

        fbkg_1 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("signal"))->getVal();
        fbkg_2 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideLo"))->getVal();
        fbkg_3 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideHi"))->getVal();
        fbkg_0 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("full"))->getVal();

        std::cout << "yields in mass windows" << std::endl;
        std::cout << sigYield.getVal()*fsig_1/fsig_0 << "\t" << bkgYield.getVal()*fbkg_1/fbkg_0 << std::endl;
        std::cout << sigYield.getVal()*(fsig_2+fsig_3)/fsig_0 << "\t" << bkgYield.getVal()*(fbkg_2+fbkg_3)/fbkg_0 << std::endl;
	
	fSigInPeak.setVal(fsig_1/fsig_0);
	fBkgInPeak.setVal(fbkg_1/fbkg_0);

	RooFitResult * result4 = data_pdf2.fitTo( *ds4Peak, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	plot(D0M, 1784., 1944., ds4, data_pdf, sig_pdfs, bkg_pdfs, "D0M_4", "m_{K#pi}");
	plot(D0LOGIPCHI2, -5., 15., *ds4Peak, data_pdf2, sig_pdfs2, bkg_pdfs2, "D0IPChi2_4", "log(IP#chi^{2})");

	print("D0params_4.dat",params);

	//Start SPlot of file4
	//SPlot only works with RooRealVar objects so need to define a new PDF
	promptYieldSP.setVal(promptYield.getVal());
	displacedYieldSP.setVal(displacedYield.getVal());
	bkgYieldSP.setVal(bkgInPeakYield.getVal());
        RooStats::SPlot* sData4 = new RooStats::SPlot("sData4","An SPlot",
                        *ds4Peak, &data_pdf2SP, RooArgList(promptYieldSP, displacedYieldSP, bkgYieldSP) );
        RooDataSet * dsSW4 = new RooDataSet("dsSW4","dsSW4",ds4Peak,*(ds4Peak->get()),0);//,"promptYieldSP_sw");

        TTree * tree_data1 = (TTree*)dsSW4->store()->tree();
        TFile * newFile1 = TFile::Open("data4_D0_sWeights.root","RECREATE");
        tree_data1->SetName("T");
        tree_data1->Write();
        newFile1->Save();
        newFile1->Close();
	//Finish SPlot

	TFile* df5 = TFile::Open("for_yandex_data_SV_5tag_new170911_DsOnly_forFit.root");
	TTree* dt5 = dynamic_cast<TTree*>(df5->Get("T"));
	RooDataSet ds5("ds5","ds5", obs, RooFit::Import(*dt5));
	RooDataSet* ds5Peak = dynamic_cast<RooDataSet*>(ds5.reduce("D0M>1844 && D0M<1884"));

	sigYield.setConstant(false);
	bkgYield.setConstant(false);
	dMass.setConstant(false);
	dWidth.setConstant(false);
	p0.setConstant(false);

	RooFitResult * result5 = data_pdf.fitTo( ds5, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	sigYield.setConstant();
	bkgYield.setConstant();
	p0.setConstant();

        fsig_1 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("signal"))->getVal();
        fsig_2 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideLo"))->getVal();
        fsig_3 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideHi"))->getVal();
        fsig_0 = sigMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("full"))->getVal();

        fbkg_1 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("signal"))->getVal();
        fbkg_2 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideLo"))->getVal();
        fbkg_3 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("sideHi"))->getVal();
        fbkg_0 = bkgMass.createIntegral(RooArgSet(D0M),RooFit::NormSet(D0M),RooFit::Range("full"))->getVal();

        std::cout << "yields in mass windows" << std::endl;
        std::cout << sigYield.getVal()*fsig_1/fsig_0 << "\t" << bkgYield.getVal()*fbkg_1/fbkg_0 << std::endl;
        std::cout << sigYield.getVal()*(fsig_2+fsig_3)/fsig_0 << "\t" << bkgYield.getVal()*(fbkg_2+fbkg_3)/fbkg_0 << std::endl;
	
	fSigInPeak.setVal(fsig_1/fsig_0);
	fBkgInPeak.setVal(fbkg_1/fbkg_0);

	RooFitResult * result6 = data_pdf2.fitTo( *ds5Peak, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	plot(D0M, 1784., 1944., ds5, data_pdf, sig_pdfs, bkg_pdfs, "D0M_5", "m_{K#pi}");
	plot(D0LOGIPCHI2, -5., 15., *ds5Peak, data_pdf2, sig_pdfs2, bkg_pdfs2, "D0IPChi2_5", "log(IP#chi^{2})");

	print("D0params_5.dat",params);

	//Start SPlot of file5
	//SPlot only works with RooRealVar objects so need to define a new PDF
	promptYieldSP.setVal(promptYield.getVal());
	displacedYieldSP.setVal(displacedYield.getVal());
	bkgYieldSP.setVal(bkgInPeakYield.getVal());
        RooStats::SPlot* sData5 = new RooStats::SPlot("sData5","An SPlot",
                        *ds5Peak, &data_pdf2SP, RooArgList(promptYieldSP, displacedYieldSP, bkgYieldSP) );
        RooDataSet * dsSW5 = new RooDataSet("dsSW5","dsSW5",ds5Peak,*(ds5Peak->get()),0);//,"promptYieldSP_sw");

        TTree * tree_data2 = (TTree*)dsSW5->store()->tree();
        TFile * newFile2 = TFile::Open("data5_D0_sWeights.root","RECREATE");
        tree_data2->SetName("T");
        tree_data2->Write();
        newFile2->Save();
        newFile2->Close();
	//Finish SPlot

}

int main() {
	fit();
}
