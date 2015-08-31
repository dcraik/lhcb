{
	TFile* f = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/gcowan/B2Kll/data/fromPatrick/kmumu_gen.root");
	TTree* t = dynamic_cast<TTree*>(f->Get("GaussDecayTupleTool/DecayTuple"));

	TFile* newfile = new TFile("/Disk/ecdf-nfs-ppe/lhcb/dcraik/kmumu_gen_addVars.root","RECREATE");
	TTree* newtree = t->CloneTree(0);

	//4-vectors for the particles
	TLorentzVector mup_4mom;
	TLorentzVector mum_4mom;
	TLorentzVector Psi_4mom;
	TLorentzVector K_4mom;
	TLorentzVector B_4mom;
	
	//4-vectors with misIDs
	TLorentzVector mum_4mom_pi2mu;
	TLorentzVector K_4mom_mu2K;

	//4-vectors for the K+mu- combination
	TLorentzVector Kmum_4mom_D;
	TLorentzVector Kmum_4mom_Psi;

	//3-vectors for boosts
	TVector3 B_boost;
	TVector3 Psi_boost;
	TVector3 mup_boost;
	TVector3 mum_boost;

	Float_t mu1Px(0.), mu1Py(0.), mu1Pz(0.), mu1Mass(0.);
	Float_t mu2Px(0.), mu2Py(0.), mu2Pz(0.), mu2Mass(0.);
	Float_t KPx(0.),   KPy(0.),   KPz(0.),   KMass(0.);

	t->SetBranchAddress("mu1Px",   &mu1Px);
	t->SetBranchAddress("mu1Py",   &mu1Py);
	t->SetBranchAddress("mu1Pz",   &mu1Pz);
	t->SetBranchAddress("mu1Mass", &mu1Mass);
	t->SetBranchAddress("mu2Px",   &mu2Px);
	t->SetBranchAddress("mu2Py",   &mu2Py);
	t->SetBranchAddress("mu2Pz",   &mu2Pz);
	t->SetBranchAddress("mu2Mass", &mu2Mass);
	t->SetBranchAddress("KPx",     &KPx);
	t->SetBranchAddress("KPy",     &KPy);
	t->SetBranchAddress("KPz",     &KPz);
	t->SetBranchAddress("KMass",   &KMass);

	Double_t costhetal(0.), qSq(0.), mKmu_D(0.), mKmu_Psi(0.);

	newtree->Branch("costhetal", &costhetal);
	newtree->Branch("qSq",       &qSq);
	newtree->Branch("mKmu_D",    &mKmu_D);
	newtree->Branch("mKmu_Psi",  &mKmu_Psi);

	Int_t n = t->GetEntries();

	Double_t piMass(0.13957018), muMass(0.1056583715);

	for(Int_t i=0; i<n; ++i) {
		if(i % 10000 == 0) std::cout << "Entry " << i << " of " << n << "..." << std::endl;
		t->GetEntry(i);

		//Get 4 momentum of each track
		K_4mom.SetXYZM(   KPx/1000,   KPy/1000,   KPz/1000,   KMass/1000);
		mup_4mom.SetXYZM( mu2Px/1000, mu2Py/1000, mu2Pz/1000, mu2Mass/1000);
		mum_4mom.SetXYZM( mu1Px/1000, mu1Py/1000, mu1Pz/1000, mu1Mass/1000);

		mum_4mom_pi2mu.SetXYZM( mu1Px/1000, mu1Py/1000, mu1Pz/1000, piMass);
		K_4mom_mu2K.SetXYZM(    KPx/1000,   KPy/1000,   KPz/1000,   muMass);

		//Make 4 momenta of composites
		Psi_4mom = mup_4mom + mum_4mom;
		B_4mom   = Psi_4mom + K_4mom;

		Kmum_4mom_D   = K_4mom      + mum_4mom_pi2mu;
		Kmum_4mom_Psi = K_4mom_mu2K + mum_4mom;

		//Get q^2 in GeV^2
		qSq = Psi_4mom.M2();

		//Keep masses in MeV for consistency
		mKmu_D   = 1000.*Kmum_4mom_D.M();
		mKmu_Psi = 1000.*Kmum_4mom_Psi.M();

		//boost into B frame
		B_boost = B_4mom.BoostVector();
		Psi_4mom.Boost(-B_boost);
		mup_4mom.Boost(-B_boost);
		
		//Boost into Psi frame
		Psi_boost = Psi_4mom.BoostVector();
		mup_4mom.Boost(-Psi_boost);
		
		//cosThetaL is the angle between the Psi and the mu+ in the Psi rest frame
		mup_boost  = mup_4mom.BoostVector();
		costhetal  = Psi_boost.Dot(mup_boost)/(Psi_boost.Mag()*mup_boost.Mag());

		newtree->Fill();

	}

	newtree->AutoSave();
	newfile->Close();
}
