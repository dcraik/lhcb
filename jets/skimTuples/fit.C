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
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include <vector>
#include <fstream>

// function to make plot of 1D projection of the fit
void plot(RooRealVar& var, double min, double max, RooDataHist& dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title) {
	TCanvas c1;
	c1.SetBottomMargin(0.19);
	RooPlot* plot = var.frame(min, max);
	dh.plotOn( plot );
	int iCol(0);
	int colours[6] = {kBlue, kRed, kGreen+2, kMagenta, kOrange, kGray};
	int fills[6] = {3245, 3454, 3644, 3205, 3495, 3690};

	RooCmdArg normCmd = RooFit::NormRange("FIT");

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
void fit(TString input_dir, TString flavour){

	// -- variables from datasets
	RooRealVar SVM(  "SVM",  "SVM",  0.,  10000.,  "MeV/#it{c}^{2}"); 
	SVM.setBins(30);
	RooRealVar SVMCor(  "SVMCor",  "SVMCor",  0.,  10000.,  "MeV/#it{c}^{2}"); 
	SVMCor.setBins(30);

	// -- variables for shape
	RooRealVar fB(  "fB",   "fB",  0.3, 0., 1.);
	RooRealVar fC(  "fC",   "fC",  0.3, 0., 1.);

	// -- yield
	RooRealVar yield(  "yield",  "", 750,    0.0,     1e5);

	// -- load histograms from simulation
	TFile* fhists = TFile::Open(input_dir+"/templates.root");
	TH1D* hSVMB = dynamic_cast<TH1D*>(fhists->Get("SVM_5"));
	TH1D* hSVMC = dynamic_cast<TH1D*>(fhists->Get("SVM_4"));
	TH1D* hSVMQ = dynamic_cast<TH1D*>(fhists->Get("SVM_0"));

	RooDataHist histSVMB("histSVMB", "histSVMB", RooArgList(SVM), hSVMB);
	RooDataHist histSVMC("histSVMC", "histSVMC", RooArgList(SVM), hSVMC);
	RooDataHist histSVMQ("histSVMQ", "histSVMQ", RooArgList(SVM), hSVMQ);

	TH1D* hSVMCorB = dynamic_cast<TH1D*>(fhists->Get("SVMCor_5"));
	TH1D* hSVMCorC = dynamic_cast<TH1D*>(fhists->Get("SVMCor_4"));
	TH1D* hSVMCorQ = dynamic_cast<TH1D*>(fhists->Get("SVMCor_0"));

	RooDataHist histSVMCorB("histSVMCorB", "histSVMCorB", RooArgList(SVMCor), hSVMCorB);
	RooDataHist histSVMCorC("histSVMCorC", "histSVMCorC", RooArgList(SVMCor), hSVMCorC);
	RooDataHist histSVMCorQ("histSVMCorQ", "histSVMCorQ", RooArgList(SVMCor), hSVMCorQ);

	// -- simulation PDFs for each category
	RooHistPdf svmB( "svmB", "svmB", RooArgSet(SVM), histSVMB ); 
	RooHistPdf svmC( "svmC", "svmC", RooArgSet(SVM), histSVMC ); 
	RooHistPdf svmQ( "svmQ", "svmQ", RooArgSet(SVM), histSVMQ ); 

	RooHistPdf svmCorB( "svmCorB", "svmCorB", RooArgSet(SVMCor), histSVMCorB ); 
	RooHistPdf svmCorC( "svmCorC", "svmCorC", RooArgSet(SVMCor), histSVMCorC ); 
	RooHistPdf svmCorQ( "svmCorQ", "svmCorQ", RooArgSet(SVMCor), histSVMCorQ ); 

	// -- N-dimensional PDFs for each component
	RooProdPdf pdfB( "pdfB",  "pdfB", RooArgList(svmB,svmCorB) );
	RooProdPdf pdfC( "pdfC",  "pdfC", RooArgList(svmC,svmCorC) );
	RooProdPdf pdfQ( "pdfQ",  "pdfQ", RooArgList(svmQ,svmCorQ) );

	// -- total PDFs (with and without normalisation)
	RooAddPdf data_pdf1( "data_pdf1",  "data_pdf1", RooArgList(pdfB, pdfC, pdfQ), RooArgList(fB, fC) );
	RooAddPdf data_pdf( "data_pdf",  "data_pdf", RooArgList(data_pdf1), RooArgList(yield) );

	// -- add all feature observables to dataset
	RooArgSet obs;
	obs.add(SVM);
	obs.add(SVMCor);

	// -- load dataset for data
	TFile* df = TFile::Open(input_dir+"/for_yandex_data_SV_"+flavour+"tag.root");
	TTree* dt = dynamic_cast<TTree*>(df->Get("T"));
	RooDataSet ds("ds","ds", obs, RooFit::Import(*dt));
	RooDataHist dh("dh","dh", obs, ds);

	// -- fit model pdf to the dataset ----------------------------------------------
	RooFitResult * result = data_pdf.fitTo( dh, RooFit::Extended(), RooFit::Save(), RooFit::NumCPU(4), RooFit::Range("FIT"));

	//make plots
	std::vector<std::string> sig_pdfs;
	sig_pdfs.push_back( "pdfB" );
	sig_pdfs.push_back( "pdfC" );
	sig_pdfs.push_back( "pdfQ" );
	std::vector<std::string> bkg_pdfs;

	plot(SVMCor, 0., 10000., dh, data_pdf, sig_pdfs, bkg_pdfs, "SVMCor_"+flavour, "SVMCor");
	plot(SVM, 0., 10000., dh, data_pdf, sig_pdfs, bkg_pdfs, "SVM_"+flavour, "SVM");

	//print parameters
	RooArgList params;
	params.add(fB);
	params.add(fC);
	print("params"+flavour+".dat",params);

}

int main() {
	fit("/tmp/dcraik/","0");
	fit("/tmp/dcraik/","4");
	fit("/tmp/dcraik/","5");
}
