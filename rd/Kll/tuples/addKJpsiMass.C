void addKJpsiMass(TString file) {
	TString fname   = "/Disk/ecdf-nfs-ppe/lhcb/dcraik/Kee/"; fname+=file; fname+=".root";
	TString outname = "/Disk/ecdf-nfs-ppe/lhcb/dcraik/Kee/"; outname+=file; outname+="_addMass.root";

	TFile* fin = TFile::Open(fname);
	TTree* tin = dynamic_cast<TTree*>(fin->Get("DecayTree"));

	TFile* fout = new TFile(outname,"RECREATE");
	TTree* tout = tin->CloneTree(0);

	Bool_t B_plus_L0ElectronDecision_TOS, B_plus_L0HadronDecision_TOS, B_plus_L0ElectronDecision_TIS, B_plus_L0HadronDecision_TIS, B_plus_L0MuonDecision_TIS, B_plus_L0PhotonDecision_TIS;
	Bool_t B_plus_Hlt1TrackAllL0Decision_TOS, B_plus_Hlt2TopoE3BodyBBDTDecision_TOS, B_plus_Hlt2TopoE2BodyBBDTDecision_TOS, B_plus_Hlt2Topo3BodyBBDTDecision_TOS, B_plus_Hlt2Topo2BodyBBDTDecision_TOS;

	Double_t K_Kst_TRACK_GhostProb, e_plus_TRACK_GhostProb, e_minus_TRACK_GhostProb, K_Kst_ProbNNk, K_Kst_PIDe, e_plus_PIDe, e_minus_PIDe;
	Double_t B_plus_M, B_plus_M02_Subst0_e2pi, B_plus_M02, J_psi_1S_M;
	Float_t B_plus_DTFM_M;
	Double_t B_plus_HOP_M, B_plus_FDCHI2_OWNPV;
	Double_t e_plus_BremMultiplicity, e_minus_BremMultiplicity;

	Float_t mKJpsi;

	//L0 trigger
	tin->SetBranchAddress("B_plus_L0ElectronDecision_TOS", &B_plus_L0ElectronDecision_TOS);
	tin->SetBranchAddress("B_plus_L0HadronDecision_TOS",   &B_plus_L0HadronDecision_TOS);
	tin->SetBranchAddress("B_plus_L0ElectronDecision_TIS", &B_plus_L0ElectronDecision_TIS);
	tin->SetBranchAddress("B_plus_L0HadronDecision_TIS",   &B_plus_L0HadronDecision_TIS);
	tin->SetBranchAddress("B_plus_L0MuonDecision_TIS",     &B_plus_L0MuonDecision_TIS);
	tin->SetBranchAddress("B_plus_L0PhotonDecision_TIS",   &B_plus_L0PhotonDecision_TIS);

	//HLT
	tin->SetBranchAddress("B_plus_Hlt1TrackAllL0Decision_TOS",     &B_plus_Hlt1TrackAllL0Decision_TOS);
	tin->SetBranchAddress("B_plus_Hlt2TopoE3BodyBBDTDecision_TOS", &B_plus_Hlt2TopoE3BodyBBDTDecision_TOS);
	tin->SetBranchAddress("B_plus_Hlt2TopoE2BodyBBDTDecision_TOS", &B_plus_Hlt2TopoE2BodyBBDTDecision_TOS);
	tin->SetBranchAddress("B_plus_Hlt2Topo3BodyBBDTDecision_TOS",  &B_plus_Hlt2Topo3BodyBBDTDecision_TOS);
	tin->SetBranchAddress("B_plus_Hlt2Topo2BodyBBDTDecision_TOS",  &B_plus_Hlt2Topo2BodyBBDTDecision_TOS);

	//ghost prob
	tin->SetBranchAddress("K_Kst_TRACK_GhostProb",                 &K_Kst_TRACK_GhostProb);
	tin->SetBranchAddress("e_plus_TRACK_GhostProb",                &e_plus_TRACK_GhostProb);
	tin->SetBranchAddress("e_minus_TRACK_GhostProb",               &e_minus_TRACK_GhostProb);
	
	//PID
	tin->SetBranchAddress("K_Kst_ProbNNk",                         &K_Kst_ProbNNk);
	tin->SetBranchAddress("K_Kst_PIDe",                            &K_Kst_PIDe);
	tin->SetBranchAddress("e_plus_PIDe",                           &e_plus_PIDe);
	tin->SetBranchAddress("e_minus_PIDe",                          &e_minus_PIDe);

	//masses
	tin->SetBranchAddress("B_plus_M",                              &B_plus_M);
	tin->SetBranchAddress("B_plus_M02_Subst0_e2pi",                &B_plus_M02_Subst0_e2pi);
	tin->SetBranchAddress("B_plus_M02",                            &B_plus_M02);
	tin->SetBranchAddress("J_psi_1S_M",                            &J_psi_1S_M);
	tin->SetBranchAddress("B_plus_DTFM_M",                         &B_plus_DTFM_M);

	//hop
	tin->SetBranchAddress("B_plus_HOP_M",                          &B_plus_HOP_M);
	tin->SetBranchAddress("B_plus_FDCHI2_OWNPV",                   &B_plus_FDCHI2_OWNPV);

	//brem
	tin->SetBranchAddress("e_plus_BremMultiplicity",               &e_plus_BremMultiplicity);
	tin->SetBranchAddress("e_minus_BremMultiplicity",              &e_minus_BremMultiplicity);

	tout->Branch("mKJpsi", &mKJpsi);

	Int_t nEntries = tin->GetEntries();

	for(Int_t i=0; i<nEntries; ++i) {
		if(i%10000==0) std::cout << "Entry " << i << " of " << nEntries << "..." << std::endl;
		tin->GetEntry(i);

		mKJpsi = B_plus_DTFM_M;

		tout->Fill();

	}

	tout->AutoSave();
	fout->Close();
}
