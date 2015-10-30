void makeTree1D(TString fname){

	TFile *newfile = new TFile(fname+".root","recreate");
	TTree *newtree = new TTree("scan","");

	newtree->ReadFile(fname+".dat","bdt/D:S:B");
	newtree->Write();
	newfile->Close();
}
