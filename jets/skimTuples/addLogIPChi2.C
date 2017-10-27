void addLogIPChi2(int flavour) {
	TString name = "for_yandex_data_SV_"; name+=flavour; name+="tag_new170911_DsOnly.root";
	TFile* fin = TFile::Open(name);
	TTree* tin = dynamic_cast<TTree*>(fin->Get("T"));

	double D0IPCHI2, DIPCHI2, DSIPCHI2, LCIPCHI2, D2K3PIIPCHI2;
	double D0LOGIPCHI2, DLOGIPCHI2, DSLOGIPCHI2, LCLOGIPCHI2, D2K3PILOGIPCHI2;

	tin->SetBranchAddress("D0IPCHI2",     &D0IPCHI2);
	tin->SetBranchAddress("DIPCHI2",      &DIPCHI2);
	tin->SetBranchAddress("DSIPCHI2",     &DSIPCHI2);
	tin->SetBranchAddress("LCIPCHI2",     &LCIPCHI2);
	tin->SetBranchAddress("D2K3PIIPCHI2", &D2K3PIIPCHI2);

	name = "for_yandex_data_SV_"; name+=flavour; name+="tag_new170911_DsOnly_forFit.root";
	TFile* fout = new TFile(name, "RECREATE");
	TTree* tout = tin->CloneTree(0);

	tout->Branch("D0LOGIPCHI2",     &D0LOGIPCHI2);
	tout->Branch("DLOGIPCHI2",      &DLOGIPCHI2);
	tout->Branch("DSLOGIPCHI2",     &DSLOGIPCHI2);
	tout->Branch("LCLOGIPCHI2",     &LCLOGIPCHI2);
	tout->Branch("D2K3PILOGIPCHI2", &D2K3PILOGIPCHI2);

	for(int i=0; i<tin->GetEntries(); ++i) {
		tin->GetEntry(i);
		D0LOGIPCHI2 = TMath::Log(D0IPCHI2);
		DLOGIPCHI2 = TMath::Log(DIPCHI2);
		DSLOGIPCHI2 = TMath::Log(DSIPCHI2);
		LCLOGIPCHI2 = TMath::Log(LCIPCHI2);
		D2K3PILOGIPCHI2 = TMath::Log(D2K3PIIPCHI2);
		tout->Fill();
	}

	tout->AutoSave();
	fout->Close();
}

void addLogIPChi2() {
	addLogIPChi2(0);
	addLogIPChi2(4);
	addLogIPChi2(5);
}
