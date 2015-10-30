Double_t getMax(Int_t n, Double_t* vals) {
	Double_t max(-999.);

	for(Int_t i=0; i<n; ++i) {
		if(vals[i] > max) max=vals[i];
	}

	return max;
}

Double_t getMin(Int_t n, Double_t* vals) {
	Double_t min(999.);

	for(Int_t i=0; i<n; ++i) {
		if(vals[i] < min) min=vals[i];
	}

	return min;
}

void makePlots() {
	gROOT->ProcessLine(".x lhcbStyle.C");
	gStyle->SetOptStat(0);

	TFile* f[19];
	TTree* t[19];
	TTree* te[19];

	for(Int_t i=0; i<19; ++i) {
		TString fname("results_"); fname+=i; fname+="_P.root";
		f[i]  = TFile::Open(fname);
		t[i]  = dynamic_cast<TTree*>(f[i]->Get("total"));
		te[i] = dynamic_cast<TTree*>(f[i]->Get("error"));
	
	}

	Double_t  xs[19] = { 0.54, 1.55, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 11.375, 12.125, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 3.55, 18.5};
	Double_t exs[19] = { 0.44, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.375,  0.375,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5, 2.45,  3.5};
	Double_t  ys[19] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t eys[19] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

	Int_t nVars = t[0]->GetListOfBranches()->GetEntries();

	for(Int_t i=0; i<nVars+2; ++i) {
            TString varName("");
	    
	    if(i==nVars) {
	        varName = "G002";

	        for(Int_t j=0; j<19; ++j) {
	        	ys[ j] = 1+2*t[j]->GetMinimum(varName);
	        	eys[j] = 2*te[j]->GetMinimum(varName);
	        }
		varName = "FH";
	    } else {
	        if(i==nVars+1) varName = "G001";
		else varName = t[0]->GetListOfBranches()->At(i)->GetName();
	        for(Int_t j=0; j<19; ++j) {
	        	ys[ j] = t[j]->GetMinimum(varName);
	        	eys[j] = te[j]->GetMinimum(varName);
	        }
	        if(i==nVars+1) varName = "AFB";
	    }

	    Double_t yMin = getMin(19,ys)-getMax(19,eys);
	    Double_t yMax = getMax(19,ys)+getMax(19,eys);

	    if(yMax-yMin == 0) continue;

	    TH2D h("h","",1,0.0,22.,1,yMin,yMax);

	    TLine l(0.0,0.0,22.0,0.0);
	    l.SetLineStyle(kDashed);
	    l.SetLineColor(kRed);
	    
	    h.GetYaxis()->SetTitle(varName);
	    if(i==nVars) h.GetYaxis()->SetTitle("F_{H}");
	    if(i==nVars+1) h.GetYaxis()->SetTitle("A_{FB}");
	    h.GetXaxis()->SetTitle("q^{2} [GeV^{2}/c^{4}]");

	    TGraphErrors g1(17,xs,   ys,   exs,   eys   );
	    TGraphErrors g2( 2,xs+17,ys+17,exs+17,eys+17);

	    g2.SetLineColor(kBlue);
	    g2.SetMarkerColor(kBlue);

	    TCanvas c;
	    h.Draw();
	    l.Draw();
	    g1.Draw("PEsame");
	    g2.Draw("PEsame");
	    c.SaveAs("plots/fromPatrick/"+varName+".pdf");
	    c.SaveAs("plots/fromPatrick/"+varName+".png");
	}

}
