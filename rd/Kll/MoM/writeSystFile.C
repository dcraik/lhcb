void writeSystFile(TString syst) {
	std::ofstream fout(syst+".dat");

	for(Int_t i=0; i<19; ++i) {
		TString fname="results_"; fname+=i; fname+="_P_"; fname+=syst; fname+=".root";
		TString fname2="results_"; fname2+=i; fname2+="_P.root";

		TFile* f = TFile::Open(fname);
		TTree* t = static_cast<TTree*>(f->Get("total"));

		TFile* f2 = TFile::Open(fname2);
		TTree* t2 = static_cast<TTree*>(f2->Get("total"));

		fout << t->GetMinimum("G001")-t2->GetMinimum("G001") << "\n" << t->GetMinimum("G002")-t2->GetMinimum("G002") << "\n" << t->GetMinimum("G003")-t2->GetMinimum("G003") << "\n" << t->GetMinimum("G004")-t2->GetMinimum("G004") << "\n" << std::endl;

		f->Close();
		f2->Close();
	}

	fout.close();
}
