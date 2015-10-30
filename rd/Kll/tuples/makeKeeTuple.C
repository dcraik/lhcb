{
	TFile* f = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/dcraik/BuKee_offline_selected_030915.root");
	TTree* t = dynamic_cast<TTree*>(f->Get("DecayTree"));

	TVector3 flight;
	TVector3 mom;

	Float_t pFlight(0.);
	Float_t ptMiss(0.);
	Float_t ptMissSq(0.);
	Float_t B_NoPsiFit_Mcorr(0.);
	Float_t B_NoPsiFit_Mhop(0.);

	Float_t B_NoPsiFit_M[100], B_ConstBNoPsiFit_J_psi_1S_M[100], B_FullFit_M[100], 
		B_NoPsiFit_PV_X[100], B_NoPsiFit_PV_Y[100], B_NoPsiFit_PV_Z[100];
	Double_t B_ENDVERTEX_X, B_ENDVERTEX_Y, B_ENDVERTEX_Z;
	Double_t B_PX, B_PY, B_PZ;

	t->SetBranchAddress("B_NoPsiFit_M",                B_NoPsiFit_M);
	t->SetBranchAddress("B_ConstBNoPsiFit_J_psi_1S_M", B_ConstBNoPsiFit_J_psi_1S_M);
	t->SetBranchAddress("B_FullFit_M",                 B_FullFit_M);

	t->SetBranchAddress("B_PX",                        &B_PX);
	t->SetBranchAddress("B_PY",                        &B_PY);
	t->SetBranchAddress("B_PZ",                        &B_PZ);

	t->SetBranchAddress("B_NoPsiFit_PV_X",             B_NoPsiFit_PV_X);
	t->SetBranchAddress("B_NoPsiFit_PV_Y",             B_NoPsiFit_PV_Y);
	t->SetBranchAddress("B_NoPsiFit_PV_Z",             B_NoPsiFit_PV_Z);

	t->SetBranchAddress("B_ENDVERTEX_X",               &B_ENDVERTEX_X);
	t->SetBranchAddress("B_ENDVERTEX_Y",               &B_ENDVERTEX_Y);
	t->SetBranchAddress("B_ENDVERTEX_Z",               &B_ENDVERTEX_Z);

	TFile* newfile = new TFile("/Disk/ecdf-nfs-ppe/lhcb/dcraik/BuKee_reduced_addMcorr.root","RECREATE");
	TTree* newtree = new TTree("DecayTree","DecayTree");

	newtree->Branch("B_NoPsiFit_M",                B_NoPsiFit_M);
	newtree->Branch("B_ConstBNoPsiFit_J_psi_1S_M", B_ConstBNoPsiFit_J_psi_1S_M);
	newtree->Branch("B_FullFit_M",                 B_FullFit_M);
	newtree->Branch("B_NoPsiFit_Mcorr",            &B_NoPsiFit_Mcorr);

	Int_t n = t->GetEntries();
	for(Int_t i=0; i<n; ++i) {
		t->GetEntry(i);

		flight.SetXYZ(B_ENDVERTEX_X-B_NoPsiFit_PV_X[0], B_ENDVERTEX_Y-B_NoPsiFit_PV_Y[0], B_ENDVERTEX_Z-B_NoPsiFit_PV_Z[0]);
		mom.SetXYZ(B_PX, B_PY, B_PZ);

		pFlight  = mom*flight/flight.Mag();
		ptMissSq = mom.Mag2() - (pFlight*pFlight);
		ptMiss   = TMath::Sqrt(ptMissSq);

		//std::cout << mom.X() << "\t" << mom.Y() << "\t" << mom.Z() << "\t" << flight.X() << "\t" << flight.Y() << "\t" << flight.Z() << std::endl;
		//std::cout << mom*flight << "\t" << flight.Mag() << std::endl;
		//std::cout << mom.Mag() << "\t" << pFlight << std::endl;

		B_NoPsiFit_Mcorr = TMath::Sqrt(B_NoPsiFit_M[0]*B_NoPsiFit_M[0] + ptMissSq) + ptMiss;

		std::cout << pFlight << "\t" << ptMissSq << "\t" << ptMiss << "\t" << B_NoPsiFit_M[0]*B_NoPsiFit_M[0] << "\t" << B_NoPsiFit_Mcorr << std::endl;

		newtree->Fill();
	}
	newtree->AutoSave();
	newfile->Close();
}
