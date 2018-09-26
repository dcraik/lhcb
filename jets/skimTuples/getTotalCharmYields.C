#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include "RooPromptShape.h"

#include <iostream>

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

bool getFittedYield(TString file, TString dType, double& yield, double& error,
	double& valPromptMean, double& valPromptWidth, double& valPromptAsym,
	double& valPromptRhoL, double& valPromptRhoR, double& valDisplacedMean, double& valDisplacedWidth, bool fix=false) {

	bool status(true);

	double massVal(0.);
	double ratioVal(3.), fracVal(0.9), alphaVal(3.), nVal(1.);
	if(dType=="D0") {
		massVal=1864.;
		ratioVal=2.85;
		fracVal=0.89;
		alphaVal=3.0;
		nVal=1.0;
	} else if(dType=="D") {
		massVal=1870.;
		ratioVal=2.90;
		fracVal=0.89;
		alphaVal=2.8;
		nVal=1.0;
	} else if(dType=="DS") {
		massVal=1968.;
		ratioVal=3.12;
		fracVal=0.90;
		alphaVal=3.0;
		nVal=1.0;
	} else if(dType=="LC") {
		massVal=2286.;
		ratioVal=3.04;
		fracVal=0.69;
		alphaVal=2.7;
		nVal=1.0;
	} else if(dType=="D2K3PI") {
		massVal=1864.;
		ratioVal=2.85;
		fracVal=0.89;
		alphaVal=3.0;
		nVal=1.0;
	}

	// -- variables from datasets
	RooRealVar DM(  dType+"M",  dType+"M",  massVal-80.,  massVal+80.,  "MeV/#it{c}^{2}"); 
	DM.setBins(30);
	RooRealVar DLOGIPCHI2(  dType+"LOGIPCHI2",  dType+"LOGIPCHI2",  -5.,  15.,  ""); 
	DLOGIPCHI2.setBins(20);
	RooRealVar DSPHIM(  "DSPHIM",  "DSPHIM",  990.,  1050.,  ""); 
	//RooRealVar JetPT(  "JetPT",  "JetPT",  50.e3,  100.e3,  ""); 

	RooRealVar dMass("dMass","mass mean",massVal,massVal-20.,massVal+20.);
	RooRealVar dWidth("dWidth","mass width",20.,5.,50.);
	//RooGaussian sigMass("sigMass","",DM,dMass,dWidth);
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

	// -- yields
	RooRealVar sigYield(  "sigYield",  "yield D",    750,    0.0,     1e5);
	RooRealVar bkgYield(  "bkgYield",  "yield comb", 750,    0.0,     1e5);

	// -- total PDF
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(sigMass,bkgMass), RooArgList(sigYield,bkgYield) );

	// -- add all feature observables to dataset
	RooRealVar DE(  dType+"E",  dType+"E",  0.,   2e6); 
	RooRealVar DPX( dType+"PX", dType+"PX", -1e5, 1e5); 
	RooRealVar DPY( dType+"PY", dType+"PY", -1e5, 1e5); 
	RooRealVar DPZ( dType+"PZ", dType+"PZ", 0.,   2e6); 

	RooArgSet obs;
	obs.add(DM);
	obs.add(DLOGIPCHI2);
	obs.add(DE);
	obs.add(DPX);
	obs.add(DPY);
	obs.add(DPZ);
	//obs.add(JetPT);
	if(dType=="DS") obs.add(DSPHIM);

	// -- load dataset for data
	TFile* df = TFile::Open("for_yandex_data_SV_"+file+"tag_testE1k_DsOnly_forFit.root");
	if(!df) return false;
	TTree* dt = dynamic_cast<TTree*>(df->Get("T"));
	if(!dt) return false;
	RooDataSet ds("ds","ds", obs, RooFit::Import(*dt));
	TString peakCut = dType; peakCut+="M>"; peakCut+=massVal-20.; peakCut+=" && "; peakCut+=dType; peakCut+="M<"; peakCut+=massVal+20.;
	TString sideCut = dType; sideCut+="M<"; sideCut+=massVal-60.; sideCut+=" || "; sideCut+=dType; sideCut+="M>"; sideCut+=massVal+60.;
	RooDataSet* dsPeak = dynamic_cast<RooDataSet*>(ds.reduce(peakCut));
	RooDataSet* dsSB   = dynamic_cast<RooDataSet*>(ds.reduce(sideCut));

	// -- fit model pdf to the dataset ----------------------------------------------
	/*RooFitResult * result =*/ data_pdf.fitTo( ds, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	//sigYield.setConstant();
	//bkgYield.setConstant();
	RooRealVar  bkgYield_mean("bkgYield_mean","bkgYield_mean", bkgYield.getVal());
	RooRealVar  bkgYield_sigma("bkgYield_sigma","bkgYield_sigma", bkgYield.getError());
	RooGaussian bkgYield_constraint("bkgYield_constraint","bkgYield_constraint", bkgYield, bkgYield_mean, bkgYield_sigma);

        DM.setRange("signal",massVal-20.,massVal+20.);
        DM.setRange("sideLo",massVal-80.,massVal-60.);
        DM.setRange("sideHi",massVal+60.,massVal+80.);
        DM.setRange("full",  massVal-80.,massVal+80.);

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
	RooRealVar promptRhoL("promptRhoL","exp L prompt",valPromptRhoL,0.5,2.);
	RooRealVar promptRhoR("promptRhoR","exp R prompt",valPromptRhoR,0.5,2.);
	RooPromptShape promptLOGIPCHI2("promptLOGIPCHI2","",DLOGIPCHI2,promptMean,promptWidth,promptAsym,promptRhoL,promptRhoR);

	RooRealVar displacedMean("displacedMean","mean displ.",valDisplacedMean,3.,12.);
	RooRealVar displacedWidth("displacedWidth","width displ.",valDisplacedWidth,1.,8.);
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

	RooFitResult * result2 = data_pdf2_full.fitTo( *dsPeak, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"), RooFit::Constrain(RooArgSet(bkgYield)));

        double fsigIPCHI2_1 = promptLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("sigIPCHI2"))->getVal();
        double fsigIPCHI2_0 = promptLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("fullIPCHI2"))->getVal();

        double f2ndIPCHI2_1 = displacedLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("sigIPCHI2"))->getVal();
        double f2ndIPCHI2_0 = displacedLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("fullIPCHI2"))->getVal();

        double fbkgIPCHI2_1 = bkgLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("sigIPCHI2"))->getVal();
        double fbkgIPCHI2_0 = bkgLOGIPCHI2.createIntegral(RooArgSet(DLOGIPCHI2),RooFit::NormSet(DLOGIPCHI2),RooFit::Range("fullIPCHI2"))->getVal();

        std::cout << "yields in IP window" << std::endl;
        std::cout << promptYield.getVal()*fsigIPCHI2_1/fsigIPCHI2_0 << "\t" << displacedYield.getVal()*f2ndIPCHI2_1/f2ndIPCHI2_0 << "\t" << bkgInPeakYield.getVal()*fbkgIPCHI2_1/fbkgIPCHI2_0 << std::endl;

	yield = promptYield.getVal();
	error = promptYield.getPropagatedError(*result2);
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

	plot(DM, massVal-80., massVal+80., ds, data_pdf, sig_pdfs, bkg_pdfs, dType+"M_"+file, "m_{K#pi}");
	plot(DLOGIPCHI2, -5., 15., *dsPeak, data_pdf2, sig_pdfs2, bkg_pdfs2, dType+"IPChi2_"+file, "log(IP#chi^{2})");

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
	print(dType+"params_"+file+".dat",params);

	////Start SPlot of file0
	////SPlot only works with RooRealVar objects so need to define a new PDF
	//RooRealVar promptYieldSP(     "promptYieldSP",    "yield prompt", promptYield.getVal(),    0.0,     1e5);
	//RooRealVar displacedYieldSP(  "displacedYieldSP", "yield displ.", displacedYield.getVal(), 0.0,     1e5);
	//RooRealVar bkgYieldSP(        "bkgYieldSP",       "yield comb",   bkgInPeakYield.getVal(), 0.0,     1e5);
	//RooAddPdf data_pdf2SP( "data_pdf2SP",  "data_pdf2SP", RooArgList(promptLOGIPCHI2,displacedLOGIPCHI2,bkgLOGIPCHI2), RooArgList(promptYieldSP,displacedYieldSP,bkgYieldSP) );
        //RooStats::SPlot* sData0 = new RooStats::SPlot("sData","An SPlot",
        //                *dsPeak, &data_pdf2SP, RooArgList(promptYieldSP, displacedYieldSP, bkgYieldSP) );
        //RooDataSet * dsSW0 = new RooDataSet("dsSW0","dsSW0",dsPeak,*(dsPeak->get()),0);//,"promptYieldSP_sw");

        //TTree * tree_data0 = (TTree*)dsSW0->store()->tree();
        //TFile * newFile0 = TFile::Open("data0_D0_sWeights.root","RECREATE");
        //tree_data0->SetName("T");
        //tree_data0->Write();
        //newFile0->Save();
        //newFile0->Close();
	////Finish SPlot

	//sigYield.setConstant(false);
	//bkgYield.setConstant(false);
	//dMass.setConstant(false);
	//dWidth.setConstant(false);
	//p0.setConstant(false);

	return status;
}

// function to get the efficiency correction for a given D type in the given file
// takes kinematics from D's in the prompt signal region and kinematic-dependent 
// efficiency from MC and PIDCalib
bool getEffCorrection(TString file, TString dType, double& totalEff, double& totalErr) {
	bool status(true);

	TFile* f = TFile::Open("for_yandex_data_SV_"+file+"tag_testE1k_DsOnly_forFit.root");
	if(!f) return false;

	TTree* t = dynamic_cast<TTree*>(f->Get("T"));
	if(!t) return false;

	double DM(0.), DLOGIPCHI2(0.), DPX(0.), DPY(0.), DPZ(0.), DE(0.);

	t->SetBranchAddress(dType+"M",         &DM        );
	t->SetBranchAddress(dType+"LOGIPCHI2", &DLOGIPCHI2);
	t->SetBranchAddress(dType+"PX",        &DPX       );
	t->SetBranchAddress(dType+"PY",        &DPY       );
	t->SetBranchAddress(dType+"PZ",        &DPZ       );
	t->SetBranchAddress(dType+"E",         &DE        );

	double massVal(0.);
	if(dType=="D0") {
		massVal=1864.;
	} else if(dType=="D") {
		dType="Dp"; // efficiency histograms use "Dp" not "D"
		massVal=1870.;
	} else if(dType=="DS") {
		dType="Ds"; // efficiency histograms use "Ds" not "DS"
		massVal=1968.;
	} else if(dType=="LC") {
		dType="Lc"; // efficiency histograms use "Lc" not "LC"
		massVal=2286.;
	} else if(dType=="D2K3PI") {
		dType="K3pi"; // efficiency histograms use "K3pi" not "D2K3PI"
		massVal=1864.;
	}

	TFile* fh = TFile::Open("efficiencies50_4.root");
	if(!fh) return false;
	TH2D* effD = dynamic_cast<TH2D*>(fh->Get("efficiency"+dType));
	TH2D* pidD = dynamic_cast<TH2D*>(fh->Get("pid"+dType+"C"));
	if(!effD || !pidD) return false;

	//std::cout << t->GetEntries() << std::endl;
	int countD(0);
	double sumInvEff(0.);
	double sumErrSq(0.);
	for(int i=0; i<t->GetEntries(); ++i) {
		t->GetEntry(i);
		if(DM<massVal-20. || DM>massVal+20. || DLOGIPCHI2>3.0) continue;
		++countD;

		TLorentzVector p4D(DPX, DPY, DPZ, DE);
		int bin1 = effD->FindBin(p4D.Pt(),p4D.Rapidity());
		int bin2 = pidD->FindBin(p4D.Pt(),p4D.Rapidity());
		double eff = effD->GetBinContent(bin1) * pidD->GetBinContent(bin2);
		double effErr = effD->GetBinError(bin2)/effD->GetBinContent(bin2);
		double pidErr = pidD->GetBinError(bin2)/pidD->GetBinContent(bin2);
		if(eff>0) {
			sumInvEff += 1/eff;
			sumErrSq += (1. + effErr*effErr + pidErr*pidErr)/(eff*eff);
		} else {
			status = false;
		}
	}
	//std::cout << countD << "\t" << sumInvEff/countD << "+/-" << TMath::Sqrt(sumErrSq)/countD << std::endl;
	totalEff = sumInvEff/countD;
	totalErr = TMath::Sqrt(sumErrSq)/countD;
	return status;
}

TString getYields(TString dType) {
	double yield0(0.), error0(0.);
	double yield4(0.), error4(0.);
	//double yield5(0.), error5(0.);
	double valPromptMean(0.), valPromptWidth(1.), valPromptAsym(0.), valPromptRhoL(1.3), valPromptRhoR(1.3), valDisplacedMean(5.), valDisplacedWidth(3.);
	double eff0(0.), err0(0.);
	double eff4(0.), err4(0.);
	//double eff5(0.), err5(0.);

	double bfD(0.), errBFD(0.);
	double ffD(0.), errFFD(0.);
	if(dType=="D0") {
		bfD = 0.0389;
		errBFD = 0.0004;
		ffD = 0.542;
		errFFD = TMath::Sqrt(.024*.024 + .007*.007);
		valPromptMean = 0.87;
		valPromptWidth = 1.09;
		valPromptAsym = -0.29;
		valPromptRhoL = 1.30;
		valPromptRhoR = 1.69;
	} else if(dType=="D") {
		bfD = 0.0898;
		errBFD = 0.0028;
		ffD = 0.225;
		errFFD = TMath::Sqrt(.010*.010 + .005*.005);
		valPromptMean = 1.04;
		valPromptWidth = 1.10;
		valPromptAsym = -0.30;
		valPromptRhoL = 1.32;
		valPromptRhoR = 1.43;
	} else if(dType=="DS") {
		bfD = 0.0227;
		errBFD = 0.0008;
		ffD = 0.092;
		errFFD = TMath::Sqrt(.008*.008 + .005*.005);
		valPromptMean = 1.03;
		valPromptWidth = 1.11;
		valPromptAsym = -0.31;
		valPromptRhoL = 1.34;
		valPromptRhoR = 1.35;
	} else if(dType=="LC") {
		bfD = 0.0635;
		errBFD = 0.0033;
		ffD = 0.057;
		errFFD = TMath::Sqrt(.006*.006 + .003*.003);
	} else if(dType=="D2K3PI") {
		bfD = 0.0811;
		errBFD = 0.0015;
		ffD = 0.542;
		errFFD = TMath::Sqrt(.024*.024 + .007*.007);
	}
	double scaleD = 1./(bfD*ffD);
	double errScaleD = scaleD*TMath::Sqrt(errBFD*errBFD/bfD/bfD + errFFD*errFFD/ffD/ffD);
	//std::cout << scaleD << "\t" << errScaleD << std::endl;
	
	if(!getFittedYield("0", dType, yield0, error0, valPromptMean, valPromptWidth, valPromptAsym, valPromptRhoL, valPromptRhoR, valDisplacedMean, valDisplacedWidth,true)) std::cout << "getFittedYield had warnings" << std::endl;
	if(!getEffCorrection("0",dType,eff0,err0)) std::cout << "getEffCorrection had warnings" << std::endl;
	
	if(!getFittedYield("4", dType, yield4, error4, valPromptMean, valPromptWidth, valPromptAsym, valPromptRhoL, valPromptRhoR, valDisplacedMean, valDisplacedWidth,true)) std::cout << "getFittedYield had warnings" << std::endl;
	if(!getEffCorrection("4",dType,eff4,err4)) std::cout << "getEffCorrection had warnings" << std::endl;
	
	//if(!getFittedYield("for_yandex_data_SV_5tag_DsOnly_forFit.root", dType, yield5, error5, valPromptMean, valPromptWidth, valPromptAsym, valPromptRhoL, valPromptRhoR, valDisplacedMean, valDisplacedWidth,true)) std::cout << "getFittedYield had warnings" << std::endl;
	//if(!getEffCorrection("for_yandex_data_SV_5tag_DsOnly_forFit.root",dType,eff5,err5)) std::cout << "getEffCorrection had warnings" << std::endl;

	//std::cout << yield0 << "+/-" << error0 << "\t" << yield4 << "+/-" << error4 << /*"\t" << yield5 << "+/-" << error5 <<*/ std::endl;
	//std::cout << eff0 << "+/-" << err0 << "\t"  << eff4 << "+/-" << err4 << /*"\t"  << eff5 << "+/-" << err5 <<*/ std::endl;

	double nCharm0(yield0*eff0*scaleD), nCharm4(yield4*eff4*scaleD);
	double errnCharm0 = nCharm0*TMath::Sqrt(err0*err0/eff0/eff0 + error0*error0/yield0/yield0 + errScaleD*errScaleD/scaleD/scaleD);
	double errnCharm4 = nCharm4*TMath::Sqrt(err4*err4/eff4/eff4 + error4*error4/yield4/yield4 + errScaleD*errScaleD/scaleD/scaleD);
	//std::cout << nCharm0 << "+/-" << errnCharm0 << "\t" << nCharm4 << "+/-" << errnCharm4 << std::endl;

	TString result;
	result.Form("%s:\t% 6.0f +/- % 5.0f\t% 6.0f +/- % 5.0f\n%.0f +/- %.0f\t%.2f +/- %.2f\t%.1f +/- %.1f\n%.0f +/- %.0f\t%.2f +/- %.2f\t%.1f +/- %.1f", dType.Data(), nCharm0, errnCharm0, nCharm4, errnCharm4, yield0, error0, eff0, err0, scaleD, errScaleD, yield4, error4, eff4, err4, scaleD, errScaleD);
	return result;

}

int main() {
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
	TString resD0=getYields("D0");
	TString resDp=getYields("D");
	TString resDs=getYields("DS");
	TString resLc=getYields("LC");
	TString resK3=getYields("D2K3PI");

	std::cout << resD0 << std::endl;
	std::cout << resDp << std::endl;
	std::cout << resDs << std::endl;
	std::cout << resLc << std::endl;
	std::cout << resK3 << std::endl;
}
