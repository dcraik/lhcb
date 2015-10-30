{
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();
	gStyle->SetOptStat(0);

	TFile* f0 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/DATA_2011_2012_LPT_PreSel_addCorrMass_NNA_TrigVeto.root");
	TFile* f1 = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/MC12_Bu2KEE_PreSel_addCorrMass_NNA_TrigVeto.root");	

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("DecayTree"));
	TTree* t1 = dynamic_cast<TTree*>(f1->Get("DecayTree"));

	TH1D* h0 = new TH1D("h0","",100,-1.,1.);
	TH1D* h1 = new TH1D("h1","",100,-1.,1.);

	t0->Draw("NN>>h0","B_DTF_PV_Mmeas>5500&&B_DTF_PV_Mmeas<5900");
	t1->Draw("NN>>h1","B_BKGCAT<30");

	h0->Scale(1./h0->Integral());
	h1->Scale(1./h1->Integral());

	h1->SetLineColor(kBlue);
	h0->SetLineColor(kRed);

	TLegend leg(0.5,0.7,0.8,0.9);
	leg.AddEntry(h1,"Signal MC","l");
	leg.AddEntry(h0,"Sideband data","l");

	h0->SetXTitle("NN output");
	h0->SetYTitle("Candidates / 0.02");

	TCanvas can;

	h0->Draw();
	h1->Draw("same");
	leg.Draw();

	can.SaveAs("NNoutput.pdf");
}
