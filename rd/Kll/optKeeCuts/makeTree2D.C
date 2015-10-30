void makeTree2D(TString fname){

	TFile *newfile = new TFile(fname+".root","recreate");
	TTree *newtree = new TTree("scan","");

	newtree->ReadFile(fname+".dat","bdt/D:pnne:S:B");
	newtree->Write();
	newfile->Close();
}
