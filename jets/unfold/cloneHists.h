#ifndef CLONEHISTS_H
#define CLONEHISTS_H

#include "TString.h"

class TH1D;
class TH2D;
class TH3D;

TH1D* cloneTH1D(TString name, TH1D* hist);
TH2D* cloneTH2D(TString name, TH2D* hist);
TH3D* cloneTH3D(TString name, TH3D* hist);

#endif
