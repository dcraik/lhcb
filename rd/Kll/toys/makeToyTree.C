void makeToyTree(TString file){

	TFile *newfile = new TFile("toys/"+file+".root","recreate");
	TTree *newtree = new TTree("toy","");

	newtree->ReadFile("toys/"+file+".dat","cosThetaL/D");
	newtree->Write();
	newfile->Close();
}
