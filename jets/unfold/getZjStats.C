{
	TFile* f0 = TFile::Open("/data/zjet/zjet_2015.root");
	TFile* f1 = TFile::Open("/data/zjet/zjet_2016.root");
	TFile* f2 = TFile::Open("/data/zjet/zjet_2017.root");
	TFile* f3 = TFile::Open("/data/zjet/zjet_2018.root");

	TTree* t0 = dynamic_cast<TTree*>(f0->Get("T"));
	TTree* t1 = dynamic_cast<TTree*>(f1->Get("T"));
	TTree* t2 = dynamic_cast<TTree*>(f2->Get("T"));
	TTree* t3 = dynamic_cast<TTree*>(f3->Get("T"));

	printf("2015 & %4d & %4d & %4d & %4d \\\\\n",
	static_cast<int>(t0->GetEntries("JetPT>15e3&&JetPT<20e3")),
	static_cast<int>(t0->GetEntries("JetPT>20e3&&JetPT<30e3")),
	static_cast<int>(t0->GetEntries("JetPT>30e3&&JetPT<50e3")),
	static_cast<int>(t0->GetEntries("JetPT>50e3")));
	printf("2016 & %4d & %4d & %4d & %4d \\\\\n",
	static_cast<int>(t1->GetEntries("JetPT>15e3&&JetPT<20e3")),
	static_cast<int>(t1->GetEntries("JetPT>20e3&&JetPT<30e3")),
	static_cast<int>(t1->GetEntries("JetPT>30e3&&JetPT<50e3")),
	static_cast<int>(t1->GetEntries("JetPT>50e3")));
	printf("2017 & %4d & %4d & %4d & %4d \\\\\n",
	static_cast<int>(t2->GetEntries("JetPT>15e3&&JetPT<20e3")),
	static_cast<int>(t2->GetEntries("JetPT>20e3&&JetPT<30e3")),
	static_cast<int>(t2->GetEntries("JetPT>30e3&&JetPT<50e3")),
	static_cast<int>(t2->GetEntries("JetPT>50e3")));
	printf("2018 & %4d & %4d & %4d & %4d \\\\\n",
	static_cast<int>(t3->GetEntries("JetPT>15e3&&JetPT<20e3")),
	static_cast<int>(t3->GetEntries("JetPT>20e3&&JetPT<30e3")),
	static_cast<int>(t3->GetEntries("JetPT>30e3&&JetPT<50e3")),
	static_cast<int>(t3->GetEntries("JetPT>50e3")));
}
