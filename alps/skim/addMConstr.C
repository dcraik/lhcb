{
	TFile* fin = TFile::Open("/data/alps/kpietapipi_addNN.root");
	TTree* tin = static_cast<TTree*>(fin->Get("T"));

	double k_e, k_px, k_py, k_pz;
	double pi_e, pi_px, pi_py, pi_pz;
	double xpip_e, xpip_px, xpip_py, xpip_pz;
	double xpim_e, xpim_px, xpim_py, xpim_pz;
	double xeta_e, xeta_px, xeta_py, xeta_pz;

	tin->SetBranchAddress("k_e",  &k_e);
	tin->SetBranchAddress("k_px", &k_px);
	tin->SetBranchAddress("k_py", &k_py);
	tin->SetBranchAddress("k_pz", &k_pz);
	tin->SetBranchAddress("pi_e",  &pi_e);
	tin->SetBranchAddress("pi_px", &pi_px);
	tin->SetBranchAddress("pi_py", &pi_py);
	tin->SetBranchAddress("pi_pz", &pi_pz);
	tin->SetBranchAddress("xpip_e",  &xpip_e);
	tin->SetBranchAddress("xpip_px", &xpip_px);
	tin->SetBranchAddress("xpip_py", &xpip_py);
	tin->SetBranchAddress("xpip_pz", &xpip_pz);
	tin->SetBranchAddress("xpim_e",  &xpim_e);
	tin->SetBranchAddress("xpim_px", &xpim_px);
	tin->SetBranchAddress("xpim_py", &xpim_py);
	tin->SetBranchAddress("xpim_pz", &xpim_pz);
	tin->SetBranchAddress("xeta_e",  &xeta_e);
	tin->SetBranchAddress("xeta_px", &xeta_px);
	tin->SetBranchAddress("xeta_py", &xeta_py);
	tin->SetBranchAddress("xeta_pz", &xeta_pz);

	TFile* fout = TFile::Open("/data/alps/kpietapipi_addNN_addMConstr.root","RECREATE");
	TTree* tout = tin->CloneTree(0);

	double b_m_cetapr(0.);
	double b_m_ceta(0.);
	double x_m_ceta(0.);

	tout->Branch("b_m_cetapr", &b_m_cetapr);
	tout->Branch("b_m_ceta",   &b_m_ceta);
	tout->Branch("x_m_ceta",   &x_m_ceta);

	TLorentzVector p4X;
	TLorentzVector p4X_etapr;
	TLorentzVector p4B;
	TLorentzVector p4B_etapr;

	TLorentzVector p4k;
	TLorentzVector p4pi;
	TLorentzVector p4Xpip;
	TLorentzVector p4Xpim;
	TLorentzVector p4Xeta;

	const double mK(493.677), mPi(139.57039), mEta(547.862), mEtapr(957.78);

	for(int i=0; i<tin->GetEntries(); ++i) {
		tin->GetEntry(i);

		p4k.SetXYZM(k_px, k_py, k_pz, mK);
		p4pi.SetXYZM(pi_px, pi_py, pi_pz, mPi);
		p4Xpip.SetXYZM(xpip_px, xpip_py, xpip_pz, mPi);
		p4Xpim.SetXYZM(xpim_px, xpim_py, xpim_pz, mPi);
		p4Xeta.SetXYZM(xeta_px, xeta_py, xeta_pz, mEta);

		p4X = p4Xpip + p4Xpim + p4Xeta;
		p4B = p4X + p4k + p4pi;

		p4X_etapr.SetXYZM(p4X.Px(),p4X.Py(),p4X.Pz(), mEtapr);
		p4B_etapr = p4X_etapr + p4k + p4pi;

		b_m_cetapr = p4B_etapr.M();
		b_m_ceta   = p4B.M();
		x_m_ceta   = p4X.M();

		tout->Fill();
	}

	tout->Write();
	fout->Close();
}
