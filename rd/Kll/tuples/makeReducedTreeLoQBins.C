void makeReducedTreeLoQBins(Int_t bin) {

	TString binStr; binStr+=bin;
	Double_t minQ(0.), maxQ(0.);

	switch(bin) {
	case 0:
		minQ = TMath::Sqrt(1.1e6);
		maxQ = TMath::Sqrt(2.e6);
		break;
	case 1:
		minQ = TMath::Sqrt(2.e6);
		maxQ = TMath::Sqrt(3.e6);
		break;
	case 2:
		minQ = TMath::Sqrt(3.e6);
		maxQ = TMath::Sqrt(4.e6);
		break;
	case 3:
		minQ = TMath::Sqrt(4.e6);
		maxQ = TMath::Sqrt(5.e6);
		break;
	case 4:
		minQ = TMath::Sqrt(5.e6);
		maxQ = TMath::Sqrt(6.e6);
		break;
	default:
		return;
	}

	TFile* f = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/gcowan/B2Kll/data/fromAlex/BuKmm.root");
	TTree* tree = dynamic_cast<TTree*>(f->Get("DecayTree"));

	TFile* fs = TFile::Open("Kmm_loQ_"+binStr+"_sWeights.root");
	TTree* stree = dynamic_cast<TTree*>(fs->Get("sWeights"));

	TFile* fn = new TFile("Kmm_loQ_"+binStr+"_reduced.root","RECREATE");
	TTree* newtree = tree->CloneTree(0);

	TLorentzVector mup_4mom;
	TLorentzVector mum_4mom;
	TLorentzVector Psi_4mom;
	TLorentzVector K_4mom;
	TLorentzVector B_4mom;
	
	TVector3 B_boost;
	TVector3 Psi_boost;
	TVector3 mup_boost;
	TVector3 mum_boost;

	Double_t K_PE(0.),   K_PX(0.),   K_PY(0.),   K_PZ(0.);
	Double_t mup_PE(0.), mup_PX(0.), mup_PY(0.), mup_PZ(0.);
	Double_t mum_PE(0.), mum_PX(0.), mum_PY(0.), mum_PZ(0.);

	Double_t B_M(0.), Psi_M(0.);

	Double_t cosThetaL(0.);

	Double_t sWeight(0.);

	newtree->Branch("cosThetaL"  ,&cosThetaL);
	newtree->Branch("sWeight"     ,&sWeight);

	tree->SetBranchAddress("B_M",  &B_M);
	tree->SetBranchAddress("Psi_M",  &Psi_M);

	tree->SetBranchAddress("Kplus_PE",   &K_PE);
	tree->SetBranchAddress("Kplus_PX",   &K_PX);
	tree->SetBranchAddress("Kplus_PY",   &K_PY);
	tree->SetBranchAddress("Kplus_PZ",   &K_PZ);
	tree->SetBranchAddress("muplus_PE",  &mup_PE);
	tree->SetBranchAddress("muplus_PX",  &mup_PX);
	tree->SetBranchAddress("muplus_PY",  &mup_PY);
	tree->SetBranchAddress("muplus_PZ",  &mup_PZ);
	tree->SetBranchAddress("muminus_PE", &mum_PE);
	tree->SetBranchAddress("muminus_PX", &mum_PX);
	tree->SetBranchAddress("muminus_PY", &mum_PY);
	tree->SetBranchAddress("muminus_PZ", &mum_PZ);



	stree->SetBranchAddress("Nsig_sw", &sWeight);

	Int_t i(0), is(0);

	Int_t n = tree->GetEntries();
	Int_t ns = stree->GetEntries();

	while( i<n && is<ns ) {
		if(is % 100 == 0) std::cout << "Entry " << is << " of " << ns << "..." << std::endl;
		do {
			++i;
			if(i==n) {//i cannot iterate past n because is would have already hit ns.
				std::cout << "This was supposed to be impossible. The two trees don't match!" << std::endl;
				return;
			}
			tree->GetEntry(i);
		} while(B_M<5150 || B_M>6000 || Psi_M<minQ || Psi_M>maxQ);

		++is;
		stree->GetEntry(is);

		//Get 4 momentum of each track
		K_4mom.SetPxPyPzE(   K_PX/1000,   K_PY/1000,   K_PZ/1000,   K_PE/1000);
		mup_4mom.SetPxPyPzE( mup_PX/1000, mup_PY/1000, mup_PZ/1000, mup_PE/1000);
		mum_4mom.SetPxPyPzE( mum_PX/1000, mum_PY/1000, mum_PZ/1000, mum_PE/1000);

		//Make 4 momenta of composites
		Psi_4mom = mup_4mom + mum_4mom;
		B_4mom   = Psi_4mom + K_4mom;

		//boost into B frame
		B_boost = B_4mom.BoostVector();
		Psi_4mom.Boost(-B_boost);
		mup_4mom.Boost(-B_boost);
		
		//Boost into Psi frame
		Psi_boost = Psi_4mom.BoostVector();
		mup_4mom.Boost(-Psi_boost);
		
		//cosThetaL is the angle between the Psi and the mu+ in the Psi rest frame
		mup_boost  = mup_4mom.BoostVector();
		cosThetaL  = Psi_boost.Dot(mup_boost)/(Psi_boost.Mag()*mup_boost.Mag());

		newtree->Fill();
	}

	newtree->AutoSave();
        fn->Close();
}
