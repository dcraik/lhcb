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

TH1D* cloneTH1D(TString name, std::vector<double> bins)
{
	TH1D* cloneHist = new TH1D(name,
			"",
			bins.size()-1,
			bins.data());
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

TH2D* cloneTH2D(TString name, TH1D* histX, TH1D* histY)
{
	if (!histX || !histY) return 0;
	TH2D* cloneHist = new TH2D(name,
			"",
			histX->GetNbinsX(),
			histX->GetXaxis()->GetXbins()->GetArray(),
			histY->GetNbinsX(),
			histY->GetXaxis()->GetXbins()->GetArray());
	return cloneHist;
}

TH2D* cloneTH2D(TString name, std::vector<double> binsX, std::vector<double> binsY)
{
	TH2D* cloneHist = new TH2D(name,
			"",
			binsX.size()-1,
			binsX.data(),
			binsY.size()-1,
			binsY.data());
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
