void makeReducedDTuple(int sample) {
	TString file1, file2, fileout;
	file1 = "for_yandex_data_SV_"; file1+=sample; file1+="tag.root";
	//file1 = "for_yandex_data_SV_"; file1+=sample; file1+="tag_MD.root";
	//file2 = "for_yandex_data_SV_"; file2+=sample; file2+="tag_MU.root";
	fileout = "for_yandex_data_SV_"; fileout+=sample; fileout+="tag_DsOnly.root";

	TChain* c = new TChain("T");
	c->AddFile(file1);
	//c->AddFile(file2);

	TFile* fout = new TFile(fileout,"RECREATE");
	TTree* tout = c->CloneTree(0);

	double D0M(0), DM(0), DSM(0), LCM(0), D2K3PIM(0);

	c->SetBranchAddress("D0M",       &D0M);
	c->SetBranchAddress("DM",        &DM);
	c->SetBranchAddress("DSM",       &DSM);
	c->SetBranchAddress("LCM",       &LCM);
	c->SetBranchAddress("D2K3PIM",   &D2K3PIM);

	int n = c->GetEntries();
	for(int i=0; i<n; ++i) {
		if(i%(n/20)==0) std::cout << i << " of " << n << std::endl;
		c->GetEntry(i);

		if(D0M>0 || DM>0 || DSM>0 || LCM>0 || D2K3PIM>0) tout->Fill();

	}

	tout->AutoSave();
	fout->Close();
}
