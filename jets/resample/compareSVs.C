{
	TFile* f1 = TFile::Open("lightjets_filtered_addVars.root");
	TTree* t1 = dynamic_cast<TTree*>(f1->Get("data"));

	TFile* f2 = TFile::Open("lightjets_filtered_addVars_resampled0.root");
	TTree* t2 = dynamic_cast<TTree*>(f2->Get("data"));

	std::vector<double>* svr_pass1;
	std::vector<double>* svr_pass2;
   	t1->SetBranchAddress("svr_pass", &svr_pass1);
   	t2->SetBranchAddress("svr_pass", &svr_pass2);

	int counters[9];

	for(int i=0; i<t1->GetEntries(); ++i) {
		t1->GetEntry(i);
		t2->GetEntry(i);

		bool in1(false), in2(false), pass1(false), pass2(false);

		if(svr_pass1->size() > 0) in1=true;
		if(svr_pass2->size() > 0) in2=true;

		for(uint i=0; i<svr_pass1->size(); ++i) {
			if(svr_pass1->at(i) == 1) pass1=true;
		}

		for(uint i=0; i<svr_pass2->size(); ++i) {
			if(svr_pass2->at(i) == 1) pass2=true;
		}
		if(!in1 && !in2)                   ++counters[0];
		if(!in1 && in2 && !pass2)          ++counters[1];
		if(!in1 && pass2)                  ++counters[2];
		if(in1 && !pass1 && !in2)          ++counters[3];
		if(in1 && !pass1 && in2 && !pass2) ++counters[4];
		if(in1 && !pass1 && pass2)         ++counters[5];
		if(pass1 && !in2)                  ++counters[6];
		if(pass1 && in2 && !pass2)         ++counters[7];
		if(pass1 && pass2)                 ++counters[8];

		if(!pass1 && pass2) std::cout << i << std::endl;
	}

	std::cout << "\t\toffline" << std::endl;
	std::cout << "online\tnoSV\tSVfail\tSVpass" << std::endl;
	for(int i=0; i<9; ++i) {
		if(i==0) std::cout << "noSV\t";
		if(i==3) std::cout << "SVfail\t";
		if(i==6) std::cout << "SVpass\t";
		std::cout << counters[i] << "\t";
		if(i%3==2) std::cout << std::endl;
	}
}
