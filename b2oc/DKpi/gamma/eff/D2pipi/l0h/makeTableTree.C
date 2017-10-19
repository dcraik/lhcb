void makeTableTree(TString file){
	TFile *newfile = new TFile("tables2012L0HadronTOS_S20realET/"+file+".root","recreate");
	TTree *newtree = new TTree("eff","");

	newtree->ReadFile("tables2012L0HadronTOS_S20realET/"+file+".dat","binStart/D:effIn:errIn:effOut:errOut");
	newtree->Write();
	newfile->Close();
}
