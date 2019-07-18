#include "cloneHists.h"

#include "TH2.h"
#include "TH3.h"

//functions to clone a histogram
TH1D* cloneTH1D(TString name, TH1D* hist)
{
	if (!hist) return 0;
	TH1D* cloneHist = new TH1D(name,
			"",
			hist->GetNbinsX(),
			hist->GetXaxis()->GetXbins()->GetArray());
	return cloneHist;
}

TH2D* cloneTH2D(TString name, TH2D* hist)
{
	if (!hist) return 0;
	TH2D* cloneHist = new TH2D(name,
			"",
			hist->GetNbinsX(),
			hist->GetXaxis()->GetXbins()->GetArray(),
			hist->GetNbinsY(),
			hist->GetYaxis()->GetXbins()->GetArray());
	return cloneHist;
}

TH3D* cloneTH3D(TString name, TH3D* hist)
{
	if (!hist) return 0;
	TH3D* cloneHist = new TH3D(name,
			"",
			hist->GetNbinsX(),
			hist->GetXaxis()->GetXbins()->GetArray(),
			hist->GetNbinsY(),
			hist->GetYaxis()->GetXbins()->GetArray(),
			hist->GetNbinsZ(),
			hist->GetZaxis()->GetXbins()->GetArray());
	return cloneHist;
}
