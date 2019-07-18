#ifndef OUTPUTFUNCTIONS_H
#define OUTPUTFUNCTIONS_H

#include <string>
#include <vector>

#include "TString.h"

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooRealVar.h"

extern TString gSaveDir;

void plotFit(RooRealVar& var, double min, double max, int nbins, RooAbsData* dh, RooAbsPdf& pdf, std::vector<std::string>& sig_pdfs, std::vector<std::string>& bkg_pdfs,  TString name, TString title);
void printParams(TString file, const RooArgList& params);

#endif
