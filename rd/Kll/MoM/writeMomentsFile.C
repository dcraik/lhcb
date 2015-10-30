{
	std::ofstream fout("moments.dat");

	for(Int_t i=0; i<19; ++i) {
		TString fname="results_"; fname+=i; fname+="_P.root";

		TFile* f = TFile::Open(fname);
		TTree* t = static_cast<TTree*>(f->Get("total"));

		fout << i << "\t" << t->GetMinimum("G000") << "\t" << t->GetMinimum("G001") << "\t" << t->GetMinimum("G002") << "\t" << t->GetMinimum("G003") << "\t" << t->GetMinimum("G004") << "\t" << t->GetMinimum("G005") << "\t" << t->GetMinimum("G006") << std::endl;

		f->Close();
	}

	fout.close();
}
