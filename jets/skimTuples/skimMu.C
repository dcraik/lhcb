{
	TFile* f = TFile::Open("/tmp/dcraik/for_yandex_data_new.root");
	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

//	double JetDijetSVDec;
//
//	t->SetBrancheAddress("JetDijetSVDec", &JetDijetSVDec);
//
//	int n = t->GetEntries();

	TFile* fout = new TFile("/tmp/dcraik/for_yandex_data_Mu.root","RECREATE");
	TTree* tout = t->CopyTree("JetDijetSVMuDec>0 || JetDijetMuMuDec>0");
	tout->Write();
	fout->Close();
}
