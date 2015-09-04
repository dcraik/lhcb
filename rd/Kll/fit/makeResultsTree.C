void makeResultsTree(TString file){

	TFile *newfile = new TFile("results_"+file+".root","recreate");
	TTree *newtree = new TTree("results","");

	newtree->ReadFile("results_"+file+".txt","i/I:a/D:b/D:NLL/d:status/I:g0/D:g1/D:g2/D:g3/D:g4/D");
	newtree->Write();
	newfile->Close();
}
