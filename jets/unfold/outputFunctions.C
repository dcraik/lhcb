#include "outputFunctions.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"

#include "RooAbsData.h"
#include "RooPlot.h"

TString gSaveDir=".";

// function to make plot of 1D projection of the fit
void plotFit(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title) {
	if(gSaveDir=="") {
		std::cout << "WARNING in plotFit: gSaveDir not set" << std::endl;
		std::cout << "                    setting to \".\"" << std::endl;
		gSaveDir=".";
	}
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

	c1.SaveAs(gSaveDir+"/"+name+"_fit.pdf");
	c1.SaveAs(gSaveDir+"/"+name+"_fit.png");
	c1.SetLogy();
	c1.SaveAs(gSaveDir+"/"+name+"_fit_log.pdf");

	delete plot;
}

// function to print parameters to a file
void printParams(TString file, const RooArgList& params) {
	if(gSaveDir=="") {
		std::cout << "WARNING in plotFit: gSaveDir not set" << std::endl;
		std::cout << "                    setting to \".\"" << std::endl;
		gSaveDir=".";
	}
	FILE * pFile = fopen((gSaveDir+"/"+file).Data(), "w");

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
