#ifndef OUTPUTFUNCTIONS_H
#define OUTPUTFUNCTIONS_H

#include <string>
#include <vector>

#include "TH1.h"
#include "TString.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooRealVar.h"

extern TString gSaveDir;

void plotFit(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf,
	     std::vector<std::string>& sig_pdfs,
	     std::vector<std::string>& bkg_pdfs,
	     TString name, TString title,
	     int typeIdx=-1,
	     RooCmdArg extraArg0=RooCmdArg::none(),
	     RooCmdArg extraArg1=RooCmdArg::none(),
	     RooCmdArg extraArg2=RooCmdArg::none());
void plotComparison(TString plotName, TString yLabel, std::vector<TH1D*> hists);

void getPrecision(double value, double error, int& precision, int& exponent);
void printParams(TString file, const RooArgList& params);
void printParamsSim(TString file, const RooArgList& params);
void printParams(TString file, RooAbsData const* data, RooAbsPdf const* pdf);
void printAllParams(TString file, RooAbsData const* data, RooAbsPdf const* pdf);
void setLHCbStyle();

#endif
