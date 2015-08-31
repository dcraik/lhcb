void makeSWeightTree(TString file){

	TFile *newfile = new TFile(file+"_sWeights.root","recreate");
	TTree *newtree = new TTree("sWeights","");

	newtree->ReadFile(file+"_sWeights.txt","B_M/D:Psi_M:Nsig_sw:Nsig_L:Nbkg_sw:Nbkg_L");
	newtree->Write();
	newfile->Close();
}
