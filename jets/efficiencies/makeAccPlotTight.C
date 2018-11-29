{
	TFile* f0 = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180903/D02Kpi_tightacc_tree.root");
	TFile* f1 = TFile::Open("/eos/user/d/dcraik/jets-tuples-new-180903/eff_09b_smooth.root");
	TFile* f2 = TFile::Open("recoEff/trkRecoEffs_fix.root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("DecayTree"));
	TH2F*  h0 = dynamic_cast<TH2F*>(f1->Get("heff"));
	TH2D*  h1 = dynamic_cast<TH2D*>(f2->Get("effPi"));
	TH2D*  h2 = dynamic_cast<TH2D*>(f2->Get("effK"));

	TH2D num("num",    "",98,2.,100.,40,2.,4.5);
	TH2D denom("denom","",98,2.,100.,40,2.,4.5);
	TH2D eff("eff",    "",98,2.,100.,40,2.,4.5);
	TH2D err("err",    "",98,2.,100.,40,2.,4.5);

	TH2D reco("reco",  "",98,2.,100.,40,2.,4.5);
	TH2D reff("reff",  "",98,2.,100.,40,2.,4.5);
	TH2D rerr("rerr",  "",98,2.,100.,40,2.,4.5);

	TH2D effcomb("effcomb",  "",98,2.,100.,40,2.,4.5);

	num.Sumw2();
	denom.Sumw2();
	reco.Sumw2();

	TCanvas c1;

	TString cuts = "K_eta_TRUE>2. && K_eta_TRUE<4.5 && K_PT_TRUE>.5 && K_P_TRUE>5. && pi_eta_TRUE>2. && pi_eta_TRUE<4.5 && pi_PT_TRUE>.5 && pi_P_TRUE>5.";

	std::cout << t0->GetEntries() << std::endl;
	std::cout << t0->GetEntries(cuts) << std::endl;

	t0->Draw("D0_y_TRUE:D0_PT_TRUE>>num",cuts);
	t0->Draw("D0_y_TRUE:D0_PT_TRUE>>denom");

	double D0_y_TRUE, D0_PT_TRUE;
	double K_eta_TRUE, K_PT_TRUE, K_P_TRUE;
	double pi_eta_TRUE, pi_PT_TRUE, pi_P_TRUE;

	t0->SetBranchAddress("D0_y_TRUE"   , &D0_y_TRUE);
	t0->SetBranchAddress("D0_PT_TRUE"  , &D0_PT_TRUE);
	t0->SetBranchAddress("K_eta_TRUE"  , &K_eta_TRUE);
	t0->SetBranchAddress("K_PT_TRUE"   , &K_PT_TRUE);
	t0->SetBranchAddress("K_P_TRUE"    , &K_P_TRUE);
	t0->SetBranchAddress("pi_eta_TRUE" , &pi_eta_TRUE);
	t0->SetBranchAddress("pi_PT_TRUE"  , &pi_PT_TRUE);
	t0->SetBranchAddress("pi_P_TRUE"   , &pi_P_TRUE);

	double recEffK(0.), recEffPi(0.);

	for(int ientry=0; ientry<t0->GetEntries(); ++ientry) {
		t0->GetEntry(ientry);
		if(!(K_eta_TRUE>2. && K_eta_TRUE<4.5 && K_PT_TRUE>.5 && K_P_TRUE>5.)) continue;
		if(!(pi_eta_TRUE>2. && pi_eta_TRUE<4.5 && pi_PT_TRUE>.5 && pi_P_TRUE>5.)) continue;

		//if(K_P_TRUE>732.) K_P_TRUE=732.;
		//if(pi_P_TRUE>732.) pi_P_TRUE=732.;
		if(K_P_TRUE>749.) K_P_TRUE=749.;
		if(pi_P_TRUE>749.) pi_P_TRUE=749.;

		//recEffK = h0->GetBinContent(h0->FindBin(K_P_TRUE,K_eta_TRUE));
		//recEffPi = h0->GetBinContent(h0->FindBin(pi_P_TRUE,pi_eta_TRUE));
		recEffK = h2->GetBinContent(h2->FindBin(K_P_TRUE*1000.,K_eta_TRUE));
		recEffPi = h1->GetBinContent(h1->FindBin(pi_P_TRUE*1000.,pi_eta_TRUE));

		reco.Fill(D0_PT_TRUE,D0_y_TRUE,recEffK*recEffPi);
	}

	eff.Divide(&num,&denom);
	reff.Divide(&reco,&num);
	for(int i=1; i<=98; ++i) {
		for(int j=1; j<=40; ++j) {
			err.SetBinContent(i,j,eff.GetBinError(i,j));
			rerr.SetBinContent(i,j,reff.GetBinError(i,j));
			effcomb.SetBinContent(i,j,eff.GetBinContent(i,j)*reff.GetBinContent(i,j));
		}
	}

	gStyle->SetOptStat(0);
	gStyle->SetPalette(1,0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs]  = { 0.00, 0.25, 0.50, 0.75, 1.00};
	Double_t reds[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 0.00};
	Double_t greens[NRGBs] = { 1.00, 0.95, 0.50, 0.00, 0.00};
	Double_t blues[NRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
	gStyle->SetNumberContours(NCont);

	eff.GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	eff.GetYaxis()->SetTitle("#it{y}(#it{D}^{0})");
	reff.GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	reff.GetYaxis()->SetTitle("#it{y}(#it{D}^{0})");
	effcomb.GetXaxis()->SetTitle("#it{p}_{T}(#it{D}^{0}) [GeV/#it{c}^{2}]");
	effcomb.GetYaxis()->SetTitle("#it{y}(#it{D}^{0})");

	eff.Draw("colz");
	c1.SaveAs("D0Acceptance.pdf");
	err.Draw("colz");
	c1.SaveAs("D0AcceptanceError.pdf");

	reff.Draw("colz");
	c1.SaveAs("D0RecoNew.pdf");
	rerr.Draw("colz");
	c1.SaveAs("D0RecoErrorNew.pdf");

	effcomb.Draw("colz");
	c1.SaveAs("D0AccRecoEff.pdf");

	TFile* fout = TFile::Open("D0AccEffNew.root","recreate");
	eff.Write();
	err.Write();
	reff.Write();
	rerr.Write();
	fout->Close();
}
