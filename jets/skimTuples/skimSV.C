{
	TChain* t = new TChain("T");
	t->Add("/tmp/dcraik/for_yandex_data_MD.root");
	t->Add("/tmp/dcraik/for_yandex_data_MU.root");
//	TFile* f = TFile::Open("/tmp/dcraik/for_yandex_data.root");
//	TTree* t = dynamic_cast<TTree*>(f->Get("T"));

//	double JetDijetSVDec;
//
//	t->SetBrancheAddress("JetDijetSVDec", &JetDijetSVDec);
//
//	int n = t->GetEntries();

	TFile* fout = new TFile("/tmp/dcraik/for_yandex_data_SV.root","RECREATE");
	TTree* tout = t->CopyTree("JetDijetSVDec>0");
	tout->Write();
	fout->Close();
}
