void makeCosThetaLPlot_P(TString bin) {
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();
	gStyle->SetOptStat(0);

	TFile* f = TFile::Open("fromPatrick/Kmm_Q"+bin+"_reduced.root");
	TTree* t = f->Get("finalTree_KMuMu");

	TH1D h("h","",50,-1.0,1.0);
	h.SetXTitle("cos #theta_{l}");
	h.SetYTitle("Candidates / 0.04");

	TCanvas c;
	t->Draw("cosThetaL>>h","sWeight","PE");
	c.SaveAs("plots/fromPatrick/cosThetaL_Q"+bin+".pdf");
}
