#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TString.h"

#include <iostream>

int main(int, char**) {

	TH1D*   hists[200];
	Double_t refs[200];
	Double_t mins[200];
	Double_t maxs[200];

	for(Int_t j=0; j<200; ++j) {
		mins[j] = 1.0;
	}
	
	TFile* f(0);
	TH1D* h(0);

	for(Int_t i=1; i<101; ++i) {
		TString fName("jackknife/"); fName+=i; fName+="/Dveto_200.root";
		f = TFile::Open(fName);
		h = dynamic_cast<TH1D*>(f->Get("efficiencyHist_18"));
		for(Int_t j=0; j<200; ++j) {
			Double_t eff = h->GetBinContent(j+1);
			if(eff>maxs[j]) maxs[j] = eff;
			if(eff<mins[j]) mins[j] = eff;
		}
		f->Close();
	}

	f = TFile::Open("Dveto_200.root");
	h = dynamic_cast<TH1D*>(f->Get("efficiencyHist_18"));

	for(Int_t j=0; j<200; ++j) {
		refs[j] = h->GetBinContent(j+1);
		TString hName("bin"); hName+=j;
		hists[j] = new TH1D(hName,"",100,mins[j],maxs[j]);
		hists[j]->SetDirectory(0);
	}

	f->Close();

	for(Int_t i=1; i<101; ++i) {
		TString fName("jackknife/"); fName+=i; fName+="/Dveto_200.root";
		f = TFile::Open(fName);
		h = dynamic_cast<TH1D*>(f->Get("efficiencyHist_18"));
		for(Int_t j=0; j<200; ++j) {
			hists[j]->Fill(h->GetBinContent(j+1));
			//std::cout << hists[j]->GetEntries() << std::endl;
//			std::cout << h->GetBinContent(j+1) << std::endl;
		}
		f->Close();
	}

	TCanvas c;
	c.Divide(5,5);

	for(Int_t i=0; i<8; ++i) {
		for(Int_t j=0; j<5; ++j) {
			for(Int_t k=0; k<5; ++k) {
				c.cd(j*5 + k + 1);
				hists[i*25 + j*5 + k]->Draw();
				Double_t x = refs[i*25 + j*5 + k];
				Double_t ymax = hists[i*25 + j*5 + k]->GetMaximum();
				TLine l(x, 0, x, ymax);
				l.SetLineColor(kRed);
				l.Draw();
			}
		}
		TString plotName("jackknifePlots_"); plotName+=i; plotName+=".pdf";
		c.SaveAs(plotName);
	}
}
