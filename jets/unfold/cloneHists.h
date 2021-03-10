#ifndef CLONEHISTS_H
#define CLONEHISTS_H

#include <vector>

#include "TString.h"

class TH1D;
class TH2D;
class TH3D;

TH1D* cloneTH1D(TString name, TH1D* hist);
TH1D* cloneTH1D(TString name, std::vector<double> bins);
TH2D* cloneTH2D(TString name, TH2D* hist);
TH2D* cloneTH2D(TString name, TH1D* histX, TH1D* histY);
TH2D* cloneTH2D(TString name, std::vector<double> binsX, std::vector<double> binsY);
TH3D* cloneTH3D(TString name, TH3D* hist);

#endif
