{
	TFile * f = TFile::Open("/Disk/ecdf-nfs-ppe/lhcb/gcowan/B2Kll/data/fromPatrick/BuKMuMu_S20_cut_2012.root");
	TTree * t = f->Get("finalTree_KMuMu");

	TFile * newfile = new TFile("fromPatrick/Kmm.root","RECREATE");
	TTree * newtree = t->CloneTree(0);

	Double_t Bplus_KMu_Jpsi(0.);
	t->SetBranchAddress("Bplus_KMu_Jpsi", &Bplus_KMu_Jpsi);

	Int_t n = t->GetEntries();

	for(Int_t i=0; i<n; ++i) {
		if(i % 50000 == 0) std::cout << "Entry " << i << " of " << n << "..." << std::endl;

		t->GetEntry(i);

		if(Bplus_KMu_Jpsi>3057&&Bplus_KMu_Jpsi<3137) continue;

		newtree->Fill();
	}
	newtree->AutoSave();
	newfile->Close();
}
