void makeTree(TString dir, Int_t i, Int_t q){

	TString fname = dir; fname+="/toy_"; fname+=i; fname+="_"; fname+=q;
	TFile *newfile = new TFile(fname+".root","recreate");
	TTree *newtree = new TTree("toy","");

	newtree->ReadFile(fname+".dat","cosThetaL/D");
	newtree->Write();
	newfile->Close();
}
