void makeToyTree(){

	TFile *newfile = new TFile("toy.root","recreate");
	TTree *newtree = new TTree("toy","");

	newtree->ReadFile("toy.txt","Mmeas/D:Mcorr");
	newtree->Write();
	newfile->Close();
}
