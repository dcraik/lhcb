{
	gStyle->SetOptStat(0000);
	gROOT->ProcessLine(".L lhcbStyle.C");
	lhcbStyle();

	TFile* f = TFile::Open("B2D0Kpi_D02KK_selBd_Dst25_PID3_NND2KK.root");
	TTree* t = ((TTree*)f->Get("DecayTree"));

	TH1D* h1 = new TH1D("h1","",50,5150,5900);
	TH1D* h2 = new TH1D("h2","",50,5150,5900);
	TH1D* h3 = new TH1D("h3","",50,5150,5900);
	TH1D* h4 = new TH1D("h4","",50,5150,5900);
	TH1D* h5 = new TH1D("h5","",50,5150,5900);
	TH1D* h6 = new TH1D("h6","",50,5150,5900);

	TCanvas c1;
	c1.Divide(3,2);
	c1.cd(1);
	t->Draw("B_D0_B_CM_M>>h1","NN>-0.7&&B_D0_B_CM_M>5150","");
	c1.cd(2);
	t->Draw("B_D0_B_CM_M>>h2","NN>-0.7&&NN<-0.4&&B_D0_B_CM_M>5150","");
	c1.cd(3);
	t->Draw("B_D0_B_CM_M>>h3","NN>-0.4&&NN<0.1&&B_D0_B_CM_M>5150","");
	c1.cd(4);
	t->Draw("B_D0_B_CM_M>>h4","NN>0.1&&NN<0.5&&B_D0_B_CM_M>5150","");
	c1.cd(5);
	t->Draw("B_D0_B_CM_M>>h5","NN>0.5&&NN<0.65&&B_D0_B_CM_M>5150","");
	c1.cd(6);
	t->Draw("B_D0_B_CM_M>>h6","NN>0.65&&B_D0_B_CM_M>5150","");

	c1.SaveAs("NNbinMassFit.pdf");
}
